subroutine 	get_RHS_EMFIE(P0,vf,ns,srcvals,zk,RHS)
implicit none

	!List of calling arguments
	integer ( kind = 4 ), intent(in) :: ns
	real ( kind = 8 ), intent(in) :: P0(3)
	complex ( kind = 8 ), intent(in) :: vf(3)
	real ( kind = 8 ), intent(in) :: srcvals(12,ns)
	complex ( kind = 8 ), intent(in) :: zk
	complex ( kind = 8 ), intent(out) :: RHS(2*ns)
	
	!List of local variables
	complex ( kind = 8 ), allocatable :: E(:,:), H(:,:)
	integer count1
	real ( kind = 8 ) ru(3),rv(3),cross_aux(3)
		allocate(E(3,ns), H(3,ns))
!	call fieldsMD(zk,P0,srcvals(1:3,:),ns,curlE,E,vf,0)
	write (*,*) 'P0',P0
	call fieldsED(zk,P0,srcvals,ns,E,H,vf,0)
	call fieldsMD(zk,P0,srcvals,ns,E,H,vf,1)
!	call point_source_vector_helmholtz(ns,P0,vf,srcvals(1:3,:),zk,E,curlE,divE)	
!	read (*,*)
	do count1=1,ns
		call orthonormalize(srcvals(4:6,count1),srcvals(10:12,count1),ru,rv)		
		RHS(count1)=-DOT_PRODUCT(rv,E(:,count1))
		RHS(ns+count1)=DOT_PRODUCT(ru,E(:,count1))
	enddo


return
end subroutine get_RHS_EMFIE

subroutine 	get_RHS_MFIE(P0,vf,ns,srcvals,zk,RHS)
implicit none

	!List of calling arguments
	integer ( kind = 4 ), intent(in) :: ns
	real ( kind = 8 ), intent(in) :: P0(3)
	complex ( kind = 8 ), intent(in) :: vf(3)
	real ( kind = 8 ), intent(in) :: srcvals(12,ns)
	complex ( kind = 8 ), intent(in) :: zk
	complex ( kind = 8 ), intent(out) :: RHS(2*ns)
	
	!List of local variables
	complex ( kind = 8 ), allocatable :: E(:,:), H(:,:)
	integer count1
	real ( kind = 8 ) ru(3),rv(3),cross_aux(3)
		allocate(E(3,ns), H(3,ns))
!	call fieldsMD(zk,P0,srcvals(1:3,:),ns,curlE,E,vf,0)
!	write (*,*) 'P0 in RHS function ',P0
	call fieldsED(zk,P0,srcvals,ns,E,H,vf,0)
	call fieldsMD(zk,P0,srcvals,ns,E,H,vf,1)
!	call point_source_vector_helmholtz(ns,P0,vf,srcvals(1:3,:),zk,E,curlE,divE)	
!	read (*,*)
	do count1=1,ns	
		call orthonormalize(srcvals(4:6,count1),srcvals(10:12,count1),ru,rv)		
		RHS(count1)=-DOT_PRODUCT(rv,H(:,count1))
		RHS(ns+count1)=DOT_PRODUCT(ru,H(:,count1))
	enddo

return
end subroutine get_RHS_MFIE


subroutine Vector_Helmholtz_targ(eps,zk,ns,source,wts,ifa_vect,a_vect,ifb_vect,&
 &b_vect,iflambda,lambda,ifrho,rho,norm_vect,ifE,E,ifcurlE,curlE,ifdivE,divE,nt,targets)
implicit none

!  This function computes the following calculation usign FMM
!
!    E=curlS_{k}[a_vect]+S_{k}[n*rho]+S_{k}[b_vect]+gradS_{k}[lambda]
!    also computes divE, curlE with appropriate flags
!
!  input:
!    eps - real *8
!      epsilon for the fmm call
!
!    zk - complex *16
!      Helmholzt parameter
!
!    ns - integer
!      number of sources
!
!    source - real *8 (3,ns)
!      location of the sources
!
!    wts - real *8 (ns)
!      weights for numerical integration
!
!    ifa_vect - integer
!      flag to consider a_vect
!
!    a_vect - *16(3,ns)
!      a vector source
!
!    ifb_vect - integer
!      flag to consider a_vect
!
!    b_vect - *16(3,ns)
!      b vector source
!
!    iflambda - integer
!      flag to consider lambda
!
!    lambda - *16(ns)
!      lambda source
!
!    ifrho - integer
!      flag to consider rho
!
!    rho - complex *16(ns)
!      rho source
!
!    norm_vect - real *8 (3,ns)
!      normal vector to the surface
!
!    ifE - integer
!      flag to compute E
!
!    ifcurlE - integer
!      flag to compute curlE
!
!    ifdivE - integer
!      flag to compute divE
!
!    targets - real *8 (3,ns)
!      location of the targets
!
!
!  output:
!
!    E - complex *16(3,nt)
!      E field
!
!    curlE - complex *16(3,nt)
!      curlE curl of E field
!
!    divE - complex *16(nt)
!      divE divergence of E
!


    !List of calling arguments
    real ( kind = 8 ), intent(in) :: eps
    complex ( kind = 8 ), intent(in) :: zk
    integer, intent(in) :: ns
	integer, intent(in) :: nt
    real ( kind = 8 ), intent(in) :: source(3,ns)
    integer, intent(in) :: ifa_vect
    complex ( kind = 8 ), intent(in) :: a_vect(3,ns)
    integer, intent(in) :: ifb_vect
    complex ( kind = 8 ), intent(in) :: b_vect(3,ns)
    integer, intent(in) :: iflambda
    complex ( kind = 8 ), intent(in) :: lambda(ns)
    integer, intent(in) :: ifrho
    complex ( kind = 8 ), intent(in) :: rho(ns)
    real ( kind = 8 ), intent(in) :: wts(ns)
    real ( kind = 8 ), intent(in) :: norm_vect(3,ns)
    integer, intent(in) :: ifE
    complex ( kind = 8 ), intent(out) :: E(3,nt)
    integer, intent(in) :: ifcurlE
    complex ( kind = 8 ), intent(out) :: curlE(3,nt)
    integer, intent(in) :: ifdivE
    complex ( kind = 8 ), intent(out) :: divE(nt)
    real ( kind = 8 ), intent(in) :: targets(3,nt)

    !List of local variables
	complex ( kind = 8 ), allocatable :: sigma_vect(:,:), dipvect_vect(:,:,:)
    complex ( kind = 8 ), allocatable :: gradE_vect(:,:,:)
    integer count1,count2,nd
    integer ier
    integer ifcharge,ifdipole,ifpot,ifgrad
	real ( kind = 8 ) pi
	
	pi=3.141592653589793238462643383d0

    !!Initialize sources
    allocate(sigma_vect(3,ns))
    allocate(dipvect_vect(3,3,ns))
    allocate(gradE_vect(3,3,nt))	

    do count1=1,ns
        sigma_vect(1,count1)=0.0d0
        sigma_vect(2,count1)=0.0d0
        sigma_vect(3,count1)=0.0d0
 
		dipvect_vect(1,1,count1)=0.0d0
		dipvect_vect(1,2,count1)=0.0d0
		dipvect_vect(1,3,count1)=0.0d0

		dipvect_vect(2,1,count1)=0.0d0
		dipvect_vect(2,2,count1)=0.0d0
		dipvect_vect(2,3,count1)=0.0d0

		dipvect_vect(3,1,count1)=0.0d0
		dipvect_vect(3,2,count1)=0.0d0
		dipvect_vect(3,3,count1)=0.0d0
    enddo

    if (ifrho.eq.1) then
        do count1=1,ns
			 sigma_vect(1,count1)=sigma_vect(1,count1)+norm_vect(1,count1)*rho(count1)
			 sigma_vect(2,count1)=sigma_vect(2,count1)+norm_vect(2,count1)*rho(count1)
			 sigma_vect(3,count1)=sigma_vect(3,count1)+norm_vect(3,count1)*rho(count1)
        enddo
    endif

    if (ifb_vect.eq.1) then
        do count1=1,ns
			 sigma_vect(1,count1)=sigma_vect(1,count1)+b_vect(1,count1)
			 sigma_vect(2,count1)=sigma_vect(2,count1)+b_vect(2,count1)
			 sigma_vect(3,count1)=sigma_vect(3,count1)+b_vect(3,count1)
        enddo
    endif

    if (iflambda.eq.1) then
        do count1=1,ns
			 dipvect_vect(1,1,count1)=dipvect_vect(1,1,count1)-lambda(count1)
			 dipvect_vect(2,2,count1)=dipvect_vect(2,2,count1)-lambda(count1)
			 dipvect_vect(3,3,count1)=dipvect_vect(3,3,count1)-lambda(count1)
        enddo
    endif

    if (ifa_vect.eq.1) then
        do count1=1,ns
!!            dipvect_x(1,count1)=dipvect_x(1,count1)+0.0d0
            dipvect_vect(1,2,count1)=dipvect_vect(1,2,count1)-a_vect(3,count1)
            dipvect_vect(1,3,count1)=dipvect_vect(1,3,count1)+a_vect(2,count1)

            dipvect_vect(2,1,count1)=dipvect_vect(2,1,count1)+a_vect(3,count1)
!!            dipvect_y(2,count1)=dipvect_y(2,count1)+0.0d0
            dipvect_vect(2,3,count1)=dipvect_vect(2,3,count1)-a_vect(1,count1)

            dipvect_vect(3,1,count1)=dipvect_vect(3,1,count1)-a_vect(2,count1)
            dipvect_vect(3,2,count1)=dipvect_vect(3,2,count1)+a_vect(1,count1)
!!            dipvect_z(3,count1)=dipvect_z(3,count1)+0.0d0

        enddo
    endif
	
	
	do count1=1,ns
		sigma_vect(1,count1)=sigma_vect(1,count1)*wts(count1)
		sigma_vect(2,count1)=sigma_vect(2,count1)*wts(count1)
		sigma_vect(3,count1)=sigma_vect(3,count1)*wts(count1)
		
		dipvect_vect(1,1,count1)=dipvect_vect(1,1,count1)*wts(count1)
		dipvect_vect(1,2,count1)=dipvect_vect(1,2,count1)*wts(count1)
		dipvect_vect(1,3,count1)=dipvect_vect(1,3,count1)*wts(count1)
		
		dipvect_vect(2,1,count1)=dipvect_vect(2,1,count1)*wts(count1)
		dipvect_vect(2,2,count1)=dipvect_vect(2,2,count1)*wts(count1)
		dipvect_vect(2,3,count1)=dipvect_vect(2,3,count1)*wts(count1)
		
		dipvect_vect(3,1,count1)=dipvect_vect(3,1,count1)*wts(count1)
		dipvect_vect(3,2,count1)=dipvect_vect(3,2,count1)*wts(count1)
		dipvect_vect(3,3,count1)=dipvect_vect(3,3,count1)*wts(count1)
	enddo

    ifcharge=1
    ifdipole=1
    ifpot=1
    ifgrad=1
!    if (ifE.eq.0) then
!        ifpot=0
!    endif
    if ((ifcurlE.eq.0).and.(ifdivE.eq.0)) then
        ifgrad=0
    endif
    if ((ifa_vect.eq.0).and.(iflambda.eq.0)) then
        ifdipole=0
    endif
    if ((ifrho.eq.0).and.(ifb_vect.eq.0)) then
        ifcharge=0
    endif
    nd=3
	if ((ifpot.eq.1).and.(ifgrad.eq.1).and.(ifcharge.eq.1).and.(ifdipole.eq.1)) then
	    call hfmm3d_t_cd_g_vec(nd,eps,zk,ns,source,sigma_vect,dipvect_vect,nt,targets,E,gradE_vect,ier)
	elseif ((ifpot.eq.1).and.(ifgrad.eq.0).and.(ifcharge.eq.1).and.(ifdipole.eq.1)) then
		call hfmm3d_t_cd_p_vec(nd,eps,zk,ns,source,sigma_vect,dipvect_vect,nt,targets,E,ier)		
	elseif ((ifpot.eq.1).and.(ifgrad.eq.1).and.(ifcharge.eq.1).and.(ifdipole.eq.0)) then
	    call hfmm3d_t_c_g_vec(nd,eps,zk,ns,source,sigma_vect,nt,targets,E,gradE_vect,ier)
	elseif ((ifpot.eq.1).and.(ifgrad.eq.0).and.(ifcharge.eq.1).and.(ifdipole.eq.0)) then
		call hfmm3d_t_c_p_vec(nd,eps,zk,ns,source,sigma_vect,nt,targets,E,ier)
	elseif ((ifpot.eq.1).and.(ifgrad.eq.1).and.(ifcharge.eq.0).and.(ifdipole.eq.1)) then
	    call hfmm3d_t_d_g_vec(nd,eps,zk,ns,source,dipvect_vect,nt,targets,E,gradE_vect,ier)
	elseif ((ifpot.eq.1).and.(ifgrad.eq.0).and.(ifcharge.eq.0).and.(ifdipole.eq.1)) then
		call hfmm3d_t_d_p_vec(nd,eps,zk,ns,source,dipvect_vect,nt,targets,E,ier)
	endif

    if (ifdivE.eq.1) then
        do count1=1,nt
			 divE(count1)=gradE_vect(1,1,count1)+gradE_vect(2,2,count1)+gradE_vect(3,3,count1)			
			 divE(count1)=divE(count1)/(4.0d0*pi)
        enddo
    endif
	if (ifE.eq.1) then
        do count1=1,nt
			 E(:,count1)=E(:,count1)/(4.0d0*pi)
        enddo
    endif

    if (ifcurlE.eq.1) then
        do count1=1,nt
		     curlE(1,count1)=gradE_vect(3,2,count1)-gradE_vect(2,3,count1)
		     curlE(2,count1)=gradE_vect(1,3,count1)-gradE_vect(3,1,count1)
		     curlE(3,count1)=gradE_vect(2,1,count1)-gradE_vect(1,2,count1)			 
			 curlE(:,count1)=curlE(:,count1)/(4.0d0*pi)
        enddo
    endif

    deallocate(sigma_vect)
    deallocate(dipvect_vect)
    deallocate(gradE_vect)


return
end subroutine Vector_Helmholtz_targ




subroutine get_thresh(srcover,ns,targs,ntarg,thresh)
implicit none

    !List of calling arguments
    integer, intent(in) :: ns,ntarg
	real ( kind = 8 ), intent(in) :: srcover(12,ns),targs(12,ntarg)
	real ( kind = 8 ), intent(out) :: thresh

    !List of local variables
	real ( kind = 8 ) xmin,xmax,ymin,ymax,zmin,zmax,boxsize,sizey,sizez
    integer i

      xmin = srcover(1,1)
      xmax = srcover(1,1) 
      ymin = srcover(2,1)
      ymax = srcover(2,1)
      zmin = srcover(3,1)
      zmax = srcover(3,1)

      do i=1,ns
        if(srcover(1,i).lt.xmin) xmin = srcover(1,i)
        if(srcover(1,i).gt.xmax) xmax = srcover(1,i)
        if(srcover(2,i).lt.ymin) ymin = srcover(2,i)
        if(srcover(2,i).gt.ymax) ymax = srcover(2,i)
        if(srcover(3,i).lt.zmin) zmin = srcover(3,i)
        if(srcover(3,i).gt.zmax) zmax = srcover(3,i)
      enddo

      do i=1,ntarg
        if(targs(1,i).lt.xmin) xmin = targs(1,i)
        if(targs(1,i).gt.xmax) xmax = targs(1,i)
        if(targs(2,i).lt.ymin) ymin = targs(2,i)
        if(targs(2,i).gt.ymax) ymax = targs(2,i)
        if(targs(3,i).lt.zmin) zmin = targs(3,i)
        if(targs(3,i).gt.zmax) zmax = targs(3,i)
      enddo
      
      boxsize = xmax-xmin
      sizey = ymax - ymin
      sizez = zmax - zmin

      if(sizey.gt.boxsize) boxsize = sizey
      if(sizez.gt.boxsize) boxsize = sizez

      thresh = 2.0d0**(-51)*boxsize


return
end subroutine



subroutine Vector_Helmholtz_targ2(eps,zk,ns,source,wts,ifa_vect,a_vect,&
 &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,norm_vect,ifE,E,ifcurlE,curlE,&
 &ifdivE,divE,nt,targets,thresh,ifdir)
implicit none

!  This function computes the following calculation usign FMM or direct
!
!    E=curlS_{k}[a_vect]+S_{k}[n*rho]+S_{k}[b_vect]+gradS_{k}[lambda]
!    also computes divE, curlE with appropriate flags
!
!  input:
!    eps - real *8
!      epsilon for the fmm call
!
!    zk - complex *16
!      Helmholzt parameter
!
!    ns - integer
!      number of sources
!
!    source - real *8 (3,ns)
!      location of the sources
!
!    wts - real *8 (ns)
!      weights for numerical integration
!
!    ifa_vect - integer
!      flag to consider a_vect
!
!    a_vect - *16(3,ns)
!      a vector source
!
!    ifb_vect - integer
!      flag to consider a_vect
!
!    b_vect - *16(3,ns)
!      b vector source
!
!    iflambda - integer
!      flag to consider lambda
!
!    lambda - *16(ns)
!      lambda source
!
!    ifrho - integer
!      flag to consider rho
!
!    rho - complex *16(ns)
!      rho source
!
!    norm_vect - real *8 (3,ns)
!      normal vector to the surface
!
!    ifE - integer
!      flag to compute E
!
!    ifcurlE - integer
!      flag to compute curlE
!
!    ifdivE - integer
!      flag to compute divE
!
!    targets - real *8 (3,ns)
!      location of the targets
!
!    thresh - real *8
!      threshold to remove the selfo interaction term
!
!    ifdir - integer
!      flag, ifdir=1 direct calculation N^2 (used to remove teh near terms)
!            ifdir=0 FMM activated
!
!  output:
!
!    E - complex *16(3,nt)
!      E field
!
!    curlE - complex *16(3,nt)
!      curlE curl of E field
!
!    divE - complex *16(nt)
!      divE divergence of E
!

  !List of calling arguments
  real ( kind = 8 ), intent(in) :: eps
  complex ( kind = 8 ), intent(in) :: zk
  integer, intent(in) :: ns
  integer, intent(in) :: nt
  real ( kind = 8 ), intent(in) :: source(3,ns)
  integer, intent(in) :: ifa_vect
  complex ( kind = 8 ), intent(in) :: a_vect(3,ns)
  integer, intent(in) :: ifb_vect
  complex ( kind = 8 ), intent(in) :: b_vect(3,ns)
  integer, intent(in) :: iflambda
  complex ( kind = 8 ), intent(in) :: lambda(ns)
  integer, intent(in) :: ifrho
  complex ( kind = 8 ), intent(in) :: rho(ns)
  real ( kind = 8 ), intent(in) :: wts(ns)
  real ( kind = 8 ), intent(in) :: norm_vect(3,ns)
  integer, intent(in) :: ifE
  complex ( kind = 8 ), intent(out) :: E(3,nt)
  integer, intent(in) :: ifcurlE
  complex ( kind = 8 ), intent(out) :: curlE(3,nt)
  integer, intent(in) :: ifdivE
  complex ( kind = 8 ), intent(out) :: divE(nt)
  real ( kind = 8 ), intent(in) :: targets(3,nt)
  real ( kind = 8 ), intent(in) :: thresh
  integer, intent(in) :: ifdir 
	

  !List of local variables
  complex ( kind = 8 ), allocatable :: sigma_vect(:,:), dipvect_vect(:,:,:)
  complex ( kind = 8 ), allocatable :: gradE_vect(:,:,:)
  integer count1,count2,nd
  integer ier
  integer ifcharge,ifdipole,ifpot,ifgrad
  real ( kind = 8 ) pi
	
  pi=3.141592653589793238462643383d0

  !!Initialize sources
  allocate(sigma_vect(3,ns))
  allocate(dipvect_vect(3,3,ns))
  allocate(gradE_vect(3,3,nt))	

  do count1=1,nt
    E(1,count1)=0.0d0
    E(2,count1)=0.0d0
    E(3,count1)=0.0d0
	gradE_vect(1,1,count1)=0.0d0
	gradE_vect(1,2,count1)=0.0d0
	gradE_vect(1,3,count1)=0.0d0
	gradE_vect(2,1,count1)=0.0d0
	gradE_vect(2,2,count1)=0.0d0
	gradE_vect(2,3,count1)=0.0d0
	gradE_vect(3,1,count1)=0.0d0
	gradE_vect(3,2,count1)=0.0d0
	gradE_vect(3,3,count1)=0.0d0
  enddo
  do count1=1,ns
    sigma_vect(1,count1)=0.0d0
    sigma_vect(2,count1)=0.0d0
    sigma_vect(3,count1)=0.0d0
 
    dipvect_vect(1,1,count1)=0.0d0
    dipvect_vect(1,2,count1)=0.0d0
    dipvect_vect(1,3,count1)=0.0d0

    dipvect_vect(2,1,count1)=0.0d0
    dipvect_vect(2,2,count1)=0.0d0
    dipvect_vect(2,3,count1)=0.0d0

    dipvect_vect(3,1,count1)=0.0d0
    dipvect_vect(3,2,count1)=0.0d0
    dipvect_vect(3,3,count1)=0.0d0
  enddo

  if (ifrho.eq.1) then
   do count1=1,ns
	sigma_vect(1,count1)=sigma_vect(1,count1)+norm_vect(1,count1)*rho(count1)
	sigma_vect(2,count1)=sigma_vect(2,count1)+norm_vect(2,count1)*rho(count1)
	sigma_vect(3,count1)=sigma_vect(3,count1)+norm_vect(3,count1)*rho(count1)
   enddo
  endif

  if (ifb_vect.eq.1) then
    do count1=1,ns
      sigma_vect(1,count1)=sigma_vect(1,count1)+b_vect(1,count1)
      sigma_vect(2,count1)=sigma_vect(2,count1)+b_vect(2,count1)
      sigma_vect(3,count1)=sigma_vect(3,count1)+b_vect(3,count1)
    enddo
  endif

  if (iflambda.eq.1) then
    do count1=1,ns
      dipvect_vect(1,1,count1)=dipvect_vect(1,1,count1)-lambda(count1)
      dipvect_vect(2,2,count1)=dipvect_vect(2,2,count1)-lambda(count1)
      dipvect_vect(3,3,count1)=dipvect_vect(3,3,count1)-lambda(count1)
    enddo
  endif

  if (ifa_vect.eq.1) then
    do count1=1,ns
!!            dipvect_x(1,count1)=dipvect_x(1,count1)+0.0d0
      dipvect_vect(1,2,count1)=dipvect_vect(1,2,count1)-a_vect(3,count1)
      dipvect_vect(1,3,count1)=dipvect_vect(1,3,count1)+a_vect(2,count1)

      dipvect_vect(2,1,count1)=dipvect_vect(2,1,count1)+a_vect(3,count1)
!!            dipvect_y(2,count1)=dipvect_y(2,count1)+0.0d0
      dipvect_vect(2,3,count1)=dipvect_vect(2,3,count1)-a_vect(1,count1)

      dipvect_vect(3,1,count1)=dipvect_vect(3,1,count1)-a_vect(2,count1)
      dipvect_vect(3,2,count1)=dipvect_vect(3,2,count1)+a_vect(1,count1)
!!            dipvect_z(3,count1)=dipvect_z(3,count1)+0.0d0

    enddo
  endif
	
  do count1=1,ns
    sigma_vect(1,count1)=sigma_vect(1,count1)*wts(count1)
    sigma_vect(2,count1)=sigma_vect(2,count1)*wts(count1)
    sigma_vect(3,count1)=sigma_vect(3,count1)*wts(count1)
		
    dipvect_vect(1,1,count1)=dipvect_vect(1,1,count1)*wts(count1)
    dipvect_vect(1,2,count1)=dipvect_vect(1,2,count1)*wts(count1)
    dipvect_vect(1,3,count1)=dipvect_vect(1,3,count1)*wts(count1)
		
    dipvect_vect(2,1,count1)=dipvect_vect(2,1,count1)*wts(count1)
    dipvect_vect(2,2,count1)=dipvect_vect(2,2,count1)*wts(count1)
    dipvect_vect(2,3,count1)=dipvect_vect(2,3,count1)*wts(count1)
		
    dipvect_vect(3,1,count1)=dipvect_vect(3,1,count1)*wts(count1)
    dipvect_vect(3,2,count1)=dipvect_vect(3,2,count1)*wts(count1)
    dipvect_vect(3,3,count1)=dipvect_vect(3,3,count1)*wts(count1)
  enddo

  ifcharge=1
  ifdipole=1
  ifpot=1
  ifgrad=1
!     if (ifE.eq.0) then
!        ifpot=0
!    endif
  if ((ifcurlE.eq.0).and.(ifdivE.eq.0)) then
    ifgrad=0
  endif
  if ((ifa_vect.eq.0).and.(iflambda.eq.0)) then
    ifdipole=0
  endif
  if ((ifrho.eq.0).and.(ifb_vect.eq.0)) then
    ifcharge=0
  endif
    
  nd=3
  if ((ifpot.eq.1).and.(ifgrad.eq.1).and.(ifcharge.eq.1).and.&
   &(ifdipole.eq.1)) then	
    if (ifdir.eq.1) then
      call h3ddirectcdg(nd,zk,source,sigma_vect,dipvect_vect,ns&
	   &,targets,nt,E,gradE_vect,thresh)
    else
      call hfmm3d_t_cd_g_vec(nd,eps,zk,ns,source,sigma_vect,dipvect_vect&
	   &,nt,targets,E,gradE_vect,ier)
    endif
  elseif ((ifpot.eq.1).and.(ifgrad.eq.0).and.(ifcharge.eq.1).and.&
   &(ifdipole.eq.1)) then
    if (ifdir.eq.1) then
      call h3ddirectcdp(nd,zk,source,sigma_vect,dipvect_vect,ns&
	   &,targets,nt,E,thresh)
    else
	  call hfmm3d_t_cd_p_vec(nd,eps,zk,ns,source,sigma_vect,dipvect_vect&
	   &,nt,targets,E,ier)
    endif
  elseif ((ifpot.eq.1).and.(ifgrad.eq.1).and.(ifcharge.eq.1).and.&
   &(ifdipole.eq.0)) then
    if (ifdir.eq.1) then
      call h3ddirectcg(nd,zk,source,sigma_vect,ns,targets,nt,E,gradE_vect&
	   &,thresh)
    else
      call hfmm3d_t_c_g_vec(nd,eps,zk,ns,source,sigma_vect,nt,targets,E,&
	   &gradE_vect,ier)
	endif	
  elseif ((ifpot.eq.1).and.(ifgrad.eq.0).and.(ifcharge.eq.1).and.&
   &(ifdipole.eq.0)) then
    if (ifdir.eq.1) then
      call h3ddirectcp(nd,zk,source,sigma_vect,ns,targets,nt,E,thresh)
    else
      call hfmm3d_t_c_p_vec(nd,eps,zk,ns,source,sigma_vect,nt,targets,E,ier)
    endif		
  elseif ((ifpot.eq.1).and.(ifgrad.eq.1).and.(ifcharge.eq.0).and.&
   &(ifdipole.eq.1)) then
    if (ifdir.eq.1) then
	  call h3ddirectdg(nd,zk,source,dipvect_vect,ns,targets,nt,E,&
	   &gradE_vect,thresh)
    else
      call hfmm3d_t_d_g_vec(nd,eps,zk,ns,source,dipvect_vect,nt,targets&
	   &,E,gradE_vect,ier)
    endif
  elseif ((ifpot.eq.1).and.(ifgrad.eq.0).and.(ifcharge.eq.0).and.&
   &(ifdipole.eq.1)) then
    if (ifdir.eq.1) then
      call h3ddirectdp(nd,zk,source,dipvect_vect,ns,targets,nt,E,thresh)
    else
      call hfmm3d_t_d_p_vec(nd,eps,zk,ns,source,dipvect_vect,nt,targets,E,ier)
    endif
  endif

  if (ifdivE.eq.1) then
    do count1=1,nt
      divE(count1)=gradE_vect(1,1,count1)+gradE_vect(2,2,count1)+&
	   &gradE_vect(3,3,count1)			
      divE(count1)=divE(count1)/(4.0d0*pi)
    enddo
  endif
  if (ifE.eq.1) then
    do count1=1,nt
      E(:,count1)=E(:,count1)/(4.0d0*pi)
    enddo
  endif

  if (ifcurlE.eq.1) then
    do count1=1,nt
      curlE(1,count1)=gradE_vect(3,2,count1)-gradE_vect(2,3,count1)
      curlE(2,count1)=gradE_vect(1,3,count1)-gradE_vect(3,1,count1)
      curlE(3,count1)=gradE_vect(2,1,count1)-gradE_vect(1,2,count1)			 
      curlE(:,count1)=curlE(:,count1)/(4.0d0*pi)
    enddo
  endif

    deallocate(sigma_vect)
    deallocate(dipvect_vect)
    deallocate(gradE_vect)


return
end subroutine Vector_Helmholtz_targ2



subroutine Vector_Helmholtz_targ2_vect(nd,eps,zk,ns,source,wts,ifa_vect,a_vect,&
 &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,norm_vect,ifE,E,ifcurlE,curlE,&
 &ifdivE,divE,nt,targets,thresh,ifdir)
implicit none

!  This function computes the following calculation usign FMM or direct
!
!    E=curlS_{k}[a_vect]+S_{k}[n*rho]+S_{k}[b_vect]+gradS_{k}[lambda]
!    also computes divE, curlE with appropriate flags
!
!  input:
!    eps - real *8
!      epsilon for the fmm call
!
!    zk - complex *16
!      Helmholzt parameter
!
!    ns - integer
!      number of source points
!
!    nd - integer
!      number of sources (number of vector sources, for Muller mainly)
!
!    source - real *8 (3,ns)
!      location of the sources
!
!    wts - real *8 (ns)
!      weights for numerical integration
!
!    ifa_vect - integer
!      flag to consider a_vect
!
!    a_vect - *16(3,ns)
!      a vector source
!
!    ifb_vect - integer
!      flag to consider a_vect
!
!    b_vect - *16(3,ns)
!      b vector source
!
!    iflambda - integer
!      flag to consider lambda
!
!    lambda - *16(ns)
!      lambda source
!
!    ifrho - integer
!      flag to consider rho
!
!    rho - complex *16(ns)
!      rho source
!
!    norm_vect - real *8 (3,ns)
!      normal vector to the surface
!
!    ifE - integer
!      flag to compute E
!
!    ifcurlE - integer
!      flag to compute curlE
!
!    ifdivE - integer
!      flag to compute divE
!
!    targets - real *8 (3,ns)
!      location of the targets
!
!    thresh - real *8
!      threshold to remove the selfo interaction term
!
!    ifdir - integer
!      flag, ifdir=1 direct calculation N^2 (used to remove teh near terms)
!            ifdir=0 FMM activated
!
!  output:
!
!    E - complex *16(3,nt)
!      E field
!
!    curlE - complex *16(3,nt)
!      curlE curl of E field
!
!    divE - complex *16(nt)
!      divE divergence of E
!

  !List of calling arguments
  real ( kind = 8 ), intent(in) :: eps
  complex ( kind = 8 ), intent(in) :: zk
  integer, intent(in) :: ns,nd
  integer, intent(in) :: nt
  real ( kind = 8 ), intent(in) :: source(3,ns)
  integer, intent(in) :: ifa_vect
  complex ( kind = 8 ), intent(in) :: a_vect(nd,3,ns)
  integer, intent(in) :: ifb_vect
  complex ( kind = 8 ), intent(in) :: b_vect(nd,3,ns)
  integer, intent(in) :: iflambda
  complex ( kind = 8 ), intent(in) :: lambda(nd,ns)
  integer, intent(in) :: ifrho
  complex ( kind = 8 ), intent(in) :: rho(nd,ns)
  real ( kind = 8 ), intent(in) :: wts(ns)
  real ( kind = 8 ), intent(in) :: norm_vect(3,ns)
  integer, intent(in) :: ifE
  complex ( kind = 8 ), intent(out) :: E(nd,3,nt)
  integer, intent(in) :: ifcurlE
  complex ( kind = 8 ), intent(out) :: curlE(nd,3,nt)
  integer, intent(in) :: ifdivE
  complex ( kind = 8 ), intent(out) :: divE(nd,nt)
  real ( kind = 8 ), intent(in) :: targets(3,nt)
  real ( kind = 8 ), intent(in) :: thresh
  integer, intent(in) :: ifdir 
	

  !List of local variables
  complex ( kind = 8 ), allocatable :: sigma_vect(:,:,:), dipvect_vect(:,:,:,:)
  complex ( kind = 8 ), allocatable :: gradE_vect(:,:,:,:)
  integer count1,count2,nd_hfmm3d
  integer ier
  integer ifcharge,ifdipole,ifpot,ifgrad
  real ( kind = 8 ) pi
	
  pi=3.141592653589793238462643383d0

  !!Initialize sources
  allocate(sigma_vect(nd,3,ns))
  allocate(dipvect_vect(nd,3,3,ns))
  allocate(gradE_vect(nd,3,3,nt))	

  do count1=1,nt
   do count2=1,nd
    E(count2,1,count1)=0.0d0
    E(count2,2,count1)=0.0d0
    E(count2,3,count1)=0.0d0
	  gradE_vect(count2,1,1,count1)=0.0d0
	  gradE_vect(count2,1,2,count1)=0.0d0
	  gradE_vect(count2,1,3,count1)=0.0d0
	  gradE_vect(count2,2,1,count1)=0.0d0
	  gradE_vect(count2,2,2,count1)=0.0d0
	  gradE_vect(count2,2,3,count1)=0.0d0
	  gradE_vect(count2,3,1,count1)=0.0d0
	  gradE_vect(count2,3,2,count1)=0.0d0
	  gradE_vect(count2,3,3,count1)=0.0d0
   enddo
  enddo
  do count1=1,ns
   do count2=1,nd
    sigma_vect(count2,1,count1)=0.0d0
    sigma_vect(count2,2,count1)=0.0d0
    sigma_vect(count2,3,count1)=0.0d0
 
    dipvect_vect(count2,1,1,count1)=0.0d0
    dipvect_vect(count2,1,2,count1)=0.0d0
    dipvect_vect(count2,1,3,count1)=0.0d0

    dipvect_vect(count2,2,1,count1)=0.0d0
    dipvect_vect(count2,2,2,count1)=0.0d0
    dipvect_vect(count2,2,3,count1)=0.0d0

    dipvect_vect(count2,3,1,count1)=0.0d0
    dipvect_vect(count2,3,2,count1)=0.0d0
    dipvect_vect(count2,3,3,count1)=0.0d0
   enddo
  enddo

  if (ifrho.eq.1) then
   do count1=1,ns
   do count2=1,nd
	sigma_vect(count2,1,count1)=sigma_vect(count2,1,count1)+norm_vect(1,count1)*rho(count2,count1)
	sigma_vect(count2,2,count1)=sigma_vect(count2,2,count1)+norm_vect(2,count1)*rho(count2,count1)
	sigma_vect(count2,3,count1)=sigma_vect(count2,3,count1)+norm_vect(3,count1)*rho(count2,count1)
   enddo
   enddo
  endif

  if (ifb_vect.eq.1) then
    do count1=1,ns
    do count2=1,nd
      sigma_vect(count2,1,count1)=sigma_vect(count2,1,count1)+b_vect(count2,1,count1)
      sigma_vect(count2,2,count1)=sigma_vect(count2,2,count1)+b_vect(count2,2,count1)
      sigma_vect(count2,3,count1)=sigma_vect(count2,3,count1)+b_vect(count2,3,count1)
    enddo
    enddo
  endif

  if (iflambda.eq.1) then
    do count1=1,ns
    do count2=1,nd
      dipvect_vect(count2,1,1,count1)=dipvect_vect(count2,1,1,count1)-lambda(count2,count1)
      dipvect_vect(count2,2,2,count1)=dipvect_vect(count2,2,2,count1)-lambda(count2,count1)
      dipvect_vect(count2,3,3,count1)=dipvect_vect(count2,3,3,count1)-lambda(count2,count1)
    enddo
    enddo
  endif

  if (ifa_vect.eq.1) then
    do count1=1,ns
    do count2=1,nd
!!            dipvect_x(1,count1)=dipvect_x(1,count1)+0.0d0
      dipvect_vect(count2,1,2,count1)=dipvect_vect(count2,1,2,count1)-a_vect(count2,3,count1)
      dipvect_vect(count2,1,3,count1)=dipvect_vect(count2,1,3,count1)+a_vect(count2,2,count1)

      dipvect_vect(count2,2,1,count1)=dipvect_vect(count2,2,1,count1)+a_vect(count2,3,count1)
!!            dipvect_y(2,count1)=dipvect_y(2,count1)+0.0d0
      dipvect_vect(count2,2,3,count1)=dipvect_vect(count2,2,3,count1)-a_vect(count2,1,count1)

      dipvect_vect(count2,3,1,count1)=dipvect_vect(count2,3,1,count1)-a_vect(count2,2,count1)
      dipvect_vect(count2,3,2,count1)=dipvect_vect(count2,3,2,count1)+a_vect(count2,1,count1)
!!            dipvect_z(3,count1)=dipvect_z(3,count1)+0.0d0
    enddo
    enddo
  endif
	
  do count1=1,ns
    do count2=1,nd
    sigma_vect(count2,1,count1)=sigma_vect(count2,1,count1)*wts(count1)
    sigma_vect(count2,2,count1)=sigma_vect(count2,2,count1)*wts(count1)
    sigma_vect(count2,3,count1)=sigma_vect(count2,3,count1)*wts(count1)
		
    dipvect_vect(count2,1,1,count1)=dipvect_vect(count2,1,1,count1)*wts(count1)
    dipvect_vect(count2,1,2,count1)=dipvect_vect(count2,1,2,count1)*wts(count1)
    dipvect_vect(count2,1,3,count1)=dipvect_vect(count2,1,3,count1)*wts(count1)
		
    dipvect_vect(count2,2,1,count1)=dipvect_vect(count2,2,1,count1)*wts(count1)
    dipvect_vect(count2,2,2,count1)=dipvect_vect(count2,2,2,count1)*wts(count1)
    dipvect_vect(count2,2,3,count1)=dipvect_vect(count2,2,3,count1)*wts(count1)
		
    dipvect_vect(count2,3,1,count1)=dipvect_vect(count2,3,1,count1)*wts(count1)
    dipvect_vect(count2,3,2,count1)=dipvect_vect(count2,3,2,count1)*wts(count1)
    dipvect_vect(count2,3,3,count1)=dipvect_vect(count2,3,3,count1)*wts(count1)
    enddo
  enddo

  ifcharge=1
  ifdipole=1
  ifpot=1
  ifgrad=1
!     if (ifE.eq.0) then
!        ifpot=0
!    endif
  if ((ifcurlE.eq.0).and.(ifdivE.eq.0)) then
    ifgrad=0
  endif
  if ((ifa_vect.eq.0).and.(iflambda.eq.0)) then
    ifdipole=0
  endif
  if ((ifrho.eq.0).and.(ifb_vect.eq.0)) then
    ifcharge=0
  endif
    
  nd_hfmm3d=3*nd
  if ((ifpot.eq.1).and.(ifgrad.eq.1).and.(ifcharge.eq.1).and.&
   &(ifdipole.eq.1)) then	
    if (ifdir.eq.1) then
      call h3ddirectcdg(nd_hfmm3d,zk,source,sigma_vect,dipvect_vect,ns&
	   &,targets,nt,E,gradE_vect,thresh)
    else
      call hfmm3d_t_cd_g_vec(nd_hfmm3d,eps,zk,ns,source,sigma_vect,dipvect_vect&
	   &,nt,targets,E,gradE_vect,ier)
    endif
  elseif ((ifpot.eq.1).and.(ifgrad.eq.0).and.(ifcharge.eq.1).and.&
   &(ifdipole.eq.1)) then
    if (ifdir.eq.1) then
      call h3ddirectcdp(nd_hfmm3d,zk,source,sigma_vect,dipvect_vect,ns&
	   &,targets,nt,E,thresh)
    else
	  call hfmm3d_t_cd_p_vec(nd_hfmm3d,eps,zk,ns,source,sigma_vect,dipvect_vect&
	   &,nt,targets,E,ier)
    endif
  elseif ((ifpot.eq.1).and.(ifgrad.eq.1).and.(ifcharge.eq.1).and.&
   &(ifdipole.eq.0)) then
    if (ifdir.eq.1) then
      call h3ddirectcg(nd_hfmm3d,zk,source,sigma_vect,ns,targets,nt,E,gradE_vect&
	   &,thresh)
    else
      call hfmm3d_t_c_g_vec(nd_hfmm3d,eps,zk,ns,source,sigma_vect,nt,targets,E,&
	   &gradE_vect,ier)
	endif	
  elseif ((ifpot.eq.1).and.(ifgrad.eq.0).and.(ifcharge.eq.1).and.&
   &(ifdipole.eq.0)) then
    if (ifdir.eq.1) then
      call h3ddirectcp(nd_hfmm3d,zk,source,sigma_vect,ns,targets,nt,E,thresh)
    else
      call hfmm3d_t_c_p_vec(nd_hfmm3d,eps,zk,ns,source,sigma_vect,nt,targets,E,ier)
    endif		
  elseif ((ifpot.eq.1).and.(ifgrad.eq.1).and.(ifcharge.eq.0).and.&
   &(ifdipole.eq.1)) then
    if (ifdir.eq.1) then
	  call h3ddirectdg(nd_hfmm3d,zk,source,dipvect_vect,ns,targets,nt,E,&
	   &gradE_vect,thresh)
    else
      call hfmm3d_t_d_g_vec(nd_hfmm3d,eps,zk,ns,source,dipvect_vect,nt,targets&
	   &,E,gradE_vect,ier)
    endif
  elseif ((ifpot.eq.1).and.(ifgrad.eq.0).and.(ifcharge.eq.0).and.&
   &(ifdipole.eq.1)) then
    if (ifdir.eq.1) then
      call h3ddirectdp(nd_hfmm3d,zk,source,dipvect_vect,ns,targets,nt,E,thresh)
    else
      call hfmm3d_t_d_p_vec(nd_hfmm3d,eps,zk,ns,source,dipvect_vect,nt,targets,E,ier)
    endif
  endif

  if (ifdivE.eq.1) then
    do count1=1,nt
    do count2=1,nd
      divE(count2,count1)=gradE_vect(count2,1,1,count1)+gradE_vect(count2,2,2,count1)+&
	   &gradE_vect(count2,3,3,count1)			
      divE(count2,count1)=divE(count2,count1)/(4.0d0*pi)
    enddo
    enddo
  endif
  if (ifE.eq.1) then
    do count1=1,nt
    do count2=1,nd
      E(count2,:,count1)=E(count2,:,count1)/(4.0d0*pi)
    enddo
    enddo
  endif

  if (ifcurlE.eq.1) then
    do count1=1,nt
    do count2=1,nd
      curlE(count2,1,count1)=gradE_vect(count2,3,2,count1)-gradE_vect(count2,2,3,count1)
      curlE(count2,2,count1)=gradE_vect(count2,1,3,count1)-gradE_vect(count2,3,1,count1)
      curlE(count2,3,count1)=gradE_vect(count2,2,1,count1)-gradE_vect(count2,1,2,count1)			 
      curlE(count2,:,count1)=curlE(count2,:,count1)/(4.0d0*pi)
    enddo
    enddo
  endif

    deallocate(sigma_vect)
    deallocate(dipvect_vect)
    deallocate(gradE_vect)


return
end subroutine Vector_Helmholtz_targ2_vect


