
subroutine Vector_Helmholtz_targ(eps,zk,ns,source,wts,ifa_vect,a_vect,ifb_vect,&
 &b_vect,iflambda,lambda,ifrho,rho,norm_vect,ifE,E,ifcurlE,curlE,ifdivE,divE,nt,targets)
implicit none

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
            dipvect_vect(2,1,count1)=dipvect_vect(2,1,count1)+a_vect(3,count1)
            dipvect_vect(1,2,count1)=dipvect_vect(1,2,count1)-a_vect(3,count1)
            dipvect_vect(1,3,count1)=dipvect_vect(1,3,count1)+a_vect(2,count1)

!!            dipvect_y(2,count1)=dipvect_y(2,count1)+0.0d0
            dipvect_vect(3,2,count1)=dipvect_vect(3,2,count1)+a_vect(1,count1)
            dipvect_vect(3,1,count1)=dipvect_vect(3,1,count1)-a_vect(2,count1)
            dipvect_vect(2,3,count1)=dipvect_vect(2,3,count1)-a_vect(1,count1)

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
        call hfmm3d_t_cd_g_vec(nd,eps,zk,ns,source,sigma_vect,dipvect_vect,nt,targets,E,gradE_vect)

    elseif ((ifpot.eq.1).and.(ifgrad.eq.0).and.(ifcharge.eq.1).and.(ifdipole.eq.1)) then
        call hfmm3d_t_cd_p_vec(nd,eps,zk,ns,source,sigma_vect,dipvect_vect,nt,targets,E)		
    elseif ((ifpot.eq.1).and.(ifgrad.eq.1).and.(ifcharge.eq.1).and.(ifdipole.eq.0)) then
        call hfmm3d_t_c_g_vec(nd,eps,zk,ns,source,sigma_vect,nt,targets,E,gradE_vect)
    elseif ((ifpot.eq.1).and.(ifgrad.eq.0).and.(ifcharge.eq.1).and.(ifdipole.eq.0)) then
        call hfmm3d_t_c_p_vec(nd,eps,zk,ns,source,sigma_vect,nt,targets,E)
    elseif ((ifpot.eq.1).and.(ifgrad.eq.1).and.(ifcharge.eq.0).and.(ifdipole.eq.1)) then
        call hfmm3d_t_d_g_vec(nd,eps,zk,ns,source,dipvect_vect,nt,targets,E,gradE_vect)
    elseif ((ifpot.eq.1).and.(ifgrad.eq.0).and.(ifcharge.eq.0).and.(ifdipole.eq.1)) then
        call hfmm3d_t_d_p_vec(nd,eps,zk,ns,source,dipvect_vect,nt,targets,E)
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



