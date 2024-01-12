


subroutine point_source_scalar_helmholtz(P0,ns,points,normals,zk,pot,dpot_dnormal)
implicit none

	!List of calling arguments
	integer ( kind = 8 ), intent(in)  :: ns
	real ( kind = 8 ), intent(in) :: P0(3)
	real ( kind = 8 ), intent(in) :: points(3,ns)
	real ( kind = 8 ), intent(in) :: normals(3,ns)
	complex ( kind = 8 ), intent(in) :: zk
	complex ( kind = 8 ), intent(out) :: pot(ns)
	complex ( kind = 8 ), intent(out) :: dpot_dnormal(ns)

    !List of local variables
	integer ( kind = 8 ) count1
	real ( kind = 8 ) r
	complex (kind = 8 ) ima
	complex ( kind = 8 ) gradpot_aux(3),aux
	
		ima=(0.0d0,1.0d0)
		
		do count1=1,ns
			r=sqrt((P0(1)-points(1,count1))**2+(P0(2)-points(2,count1))**2+(P0(3)-points(3,count1))**2)
			pot(count1)=exp(ima*zk*r)/r
			aux=(-1.0d0+ima*zk*r)/r**2*pot(count1)
			gradpot_aux(1)=aux*(points(1,count1)-P0(1))
			gradpot_aux(2)=aux*(points(2,count1)-P0(2))
			gradpot_aux(3)=aux*(points(3,count1)-P0(3))
			dpot_dnormal(count1)=gradpot_aux(1)*normals(1,count1)+gradpot_aux(2)*normals(2,count1)+gradpot_aux(3)*normals(3,count1)
		enddo

return
end subroutine point_source_scalar_helmholtz


subroutine point_source_scalar_helmholtz2(P0,qv,ns,points,zk,pot,gradpot,init)
implicit none

	!List of calling arguments
	integer ( kind = 8 ), intent(in)  :: ns,init
	real ( kind = 8 ), intent(in) :: P0(3)
	real ( kind = 8 ), intent(in) :: points(3,ns)
	complex ( kind = 8 ), intent(in) :: zk,qv
	complex ( kind = 8 ), intent(out) :: pot(ns)
	complex ( kind = 8 ), intent(out) :: gradpot(3,ns)

    !List of local variables
	integer ( kind = 8 ) count1
	real ( kind = 8 ) r
	complex (kind = 8 ) ima
	complex ( kind = 8 ) aux
	
		ima=(0.0d0,1.0d0)
		
		do count1=1,ns
		  if (init.eq.0) then
			pot(count1)=0.0d0
			gradpot(1,count1)=0.0d0
			gradpot(2,count1)=0.0d0
			gradpot(3,count1)=0.0d0
		  endif
			r=sqrt((P0(1)-points(1,count1))**2+(P0(2)-points(2,count1))**2+(P0(3)-points(3,count1))**2)
			pot(count1)=pot(count1)+qv*exp(ima*zk*r)/r
			aux=(-1.0d0+ima*zk*r)/r**3*exp(ima*zk*r)
			gradpot(1,count1)=gradpot(1,count1)+qv*aux*(points(1,count1)-P0(1))
			gradpot(2,count1)=gradpot(2,count1)+qv*aux*(points(2,count1)-P0(2))
			gradpot(3,count1)=gradpot(3,count1)+qv*aux*(points(3,count1)-P0(3))
		enddo

return
end subroutine point_source_scalar_helmholtz2

subroutine point_source_vector_helmholtz(ns,P0,vf,points,zk,E,curlE,divE)
implicit none
! a small electric dipole; that is: E=vf*exp(ima*zk*r)/r 

	!List of calling arguments
	integer ( kind = 8 ), intent(in) :: ns
	real ( kind = 8 ), intent(in) :: P0(3)
	complex ( kind = 8 ), intent(in) :: vf(3)
	real ( kind = 8 ), intent(in) :: points(3,ns)
	complex ( kind = 8 ), intent(in) :: zk
	complex ( kind = 8 ), intent(out) :: E(3,ns)
	complex ( kind = 8 ), intent(out) :: curlE(3,ns)
	complex ( kind = 8 ), intent(out) :: divE(ns)
	
    !List of local variables
	integer ( kind = 8 ) count1
	real ( kind = 8 ) dx,dy,dz,r
	complex (kind = 8 ) ima
	complex (kind = 8 ) R1
	
		ima=(0.0d0,1.0d0)
		do count1=1,ns
			dx=points(1,count1)-P0(1)
			dy=points(2,count1)-P0(2)
			dz=points(3,count1)-P0(3)
			r=dsqrt(dx**2+dy**2+dz**2)
			R1=(ima*zk/r**2-1.0d0/r**3)
			E(1,count1)=vf(1)*exp(ima*zk*r)/r
			E(2,count1)=vf(2)*exp(ima*zk*r)/r
			E(3,count1)=vf(3)*exp(ima*zk*r)/r
			divE(count1)=(vf(1)*dx+vf(2)*dy+vf(3)*dz)*R1*exp(ima*zk*r)
			curlE(1,count1)=(dy*vf(3)-dz*vf(2))*R1*exp(ima*zk*r)
			curlE(2,count1)=(dz*vf(1)-dx*vf(3))*R1*exp(ima*zk*r)
			curlE(3,count1)=(dx*vf(2)-dy*vf(1))*R1*exp(ima*zk*r)
		enddo
return
end subroutine point_source_vector_helmholtz






subroutine fieldsED(zk, P0, points, n, E, H, vf, init)
  implicit none
  !
  ! Generate the E&M field due to an electric dipole, i.e. the fields
  ! are generated as
  !
  !     H = curl(A)
  !     E = 1/(-ima*zk) * curl(curl(A))
  !
  ! where  A = vf*exp(ima*zk*r)/r
	
  ! List of calling arguments
  integer *8, intent(in) :: n,init
  real ( kind = 8 ), intent(in) :: P0(3),points(12,n)
  complex ( kind = 8 ), intent(in) :: vf(3),zk
  complex ( kind = 8 ), intent(out) :: E(3,n),H(3,n)
	
  ! List of local variables
  real ( kind = 8 ) dx,dy,dz,r
  complex ( kind = 8 ) R1,R2,au1,au2,ima
  integer *8 i
	
  ima = (0.0d0,1.0d0)
  do i=1,n
     if (init .eq. 0) then
        E(1,i)=0
        E(2,i)=0
        E(3,i)=0
        H(1,i)=0
        H(2,i)=0
        H(3,i)=0
     endif

     dx=points(1,i)-P0(1)
     dy=points(2,i)-P0(2)
     dz=points(3,i)-P0(3)
     r=sqrt(dx**2+dy**2+dz**2)

     R1=(ima*zk/r**2-1/r**3)
     R2=((ima*zk)**2/r**3-3*ima*zk/r**4+3/r**5)

     au1=(zk**2/r+R1)*exp(ima*zk*r)/(-ima*zk)
     au2=dx*vf(1)+dy*vf(2)
     au2=au2+dz*vf(3)
     au2=au2*R2*exp(ima*zk*r)/(-ima*zk)

     E(1,i)=E(1,i)+(vf(1)*au1+dx*au2)
     E(2,i)=E(2,i)+(vf(2)*au1+dy*au2)
     E(3,i)=E(3,i)+(vf(3)*au1+dz*au2)

     H(1,i)=H(1,i)+(dy*vf(3)-dz*vf(2))*R1*exp(ima*zk*r)
     H(2,i)=H(2,i)+(dz*vf(1)-dx*vf(3))*R1*exp(ima*zk*r)
     H(3,i)=H(3,i)+(dx*vf(2)-dy*vf(1))*R1*exp(ima*zk*r)
  enddo

  return
end subroutine fieldsED





subroutine fieldsMD(zk,P0,points,n,E,H,vf,initial)
implicit none
! a small magnetic dipole; that is: F=vf*exp(ima*zk*r)/r
! E=curlF
! H=1/(ima*zk)*curlcurlF

	!List of calling arguments
	integer *8, intent(in) :: n,initial
	real ( kind = 8 ), intent(in) :: P0(3),points(12,n)
	complex ( kind = 8 ), intent(in) :: vf(3),zk
	complex ( kind = 8 ), intent(out) :: E(3,n),H(3,n)

	!List of local variables
	real ( kind = 8 ) dx,dy,dz,r
	complex ( kind = 8 ) R1,R2,au1,au2,ima
	integer *8 i

	ima = (0.0d0,1.0d0)

	do i=1,n
		if (initial .eq. 0) then
			E(1,i)=0
			E(2,i)=0
			E(3,i)=0
			H(1,i)=0
			H(2,i)=0
			H(3,i)=0
       endif
		dx=points(1,i)-P0(1)
		dy=points(2,i)-P0(2)
		dz=points(3,i)-P0(3)
		r=sqrt(dx**2+dy**2+dz**2)
		R1=(ima*zk/r**2-1/r**3)
		R2=((ima*zk)**2/r**3-3*ima*zk/r**4+3/r**5)
		au1=(zk**2/r+R1)*exp(ima*zk*r)/(ima*zk)
		au2=dx*vf(1)+dy*vf(2)
		au2=au2+dz*vf(3)
		au2=au2*R2*exp(ima*zk*r)/(ima*zk)
		H(1,i)=H(1,i)+(vf(1)*au1+dx*au2)
		H(2,i)=H(2,i)+(vf(2)*au1+dy*au2)
		H(3,i)=H(3,i)+(vf(3)*au1+dz*au2)
		E(1,i)=E(1,i)+(dy*vf(3)-dz*vf(2))*R1*exp(ima*zk*r)
		E(2,i)=E(2,i)+(dz*vf(1)-dx*vf(3))*R1*exp(ima*zk*r)
		E(3,i)=E(3,i)+(dx*vf(2)-dy*vf(1))*R1*exp(ima*zk*r)
	enddo

return
end subroutine fieldsMD





subroutine fieldsPW(zk,xyz,n,nsource,E,H)
! a plane wave; that is: E=exp(ima*zk*z) ux
! H=exp(ima*zk*z) uy

	integer *8, intent(in) :: n
	real ( kind = 8 ), intent(in) :: xyz(3,n)
	complex ( kind = 8 ), intent(in) :: zk
	complex *16 E(3,n),H(3,n)
	complex *16 ima
	integer *8 i
	data ima/(0.0d0,1.0d0)/

      do i=1,n
        H(1,i)=0
        H(2,i)=exp(ima*zk*xyz(3,i))
        H(3,i)=0
        E(1,i)=exp(ima*zk*xyz(3,i))
        E(2,i)=0
        E(3,i)=0
      enddo

return
end subroutine fieldsPW





subroutine ScalarPW(zk,xyz,n,pot)
! a plane wave; that is: E=exp(ima*zk*z) ux
! H=exp(ima*zk*z) uy

	integer *8, intent(in) :: n
	real ( kind = 8 ), intent(in) :: xyz(3,n)
	complex ( kind = 8 ), intent(in) :: zk
	complex *16 pot(n)
	complex *16 ima
	real ( kind = 8 ) phi,pheta,kx,ky,kz
	integer *8 i
	data ima/(0.0d0,1.0d0)/
	real ( kind = 8 ) pi
	
	pi=3.1415926535897932384626433832795028841971d0


		!Two angles for the arrival direction
		phi=190  
		theta=70
		
		phi=phi/180*pi
		theta=theta/180*pi
		
		
		kz=cos(theta)
		ky=sin(theta)*sin(phi)
		kx=sin(theta)*cos(phi)

!      do i=1,n
!        pot(i)=exp(ima*zk*xyz(3,i))
!      enddo

!      do i=1,n
!        pot(i)=exp(ima*zk*(-.707d0*xyz(3,i)+.707d0*xyz(1,i)))
!      enddo

	  do i=1,n
	    pot(i)=exp(-ima*zk*(kx*xyz(1,i)+ky*xyz(2,i)+kz*xyz(3,i)))
	  enddo

return
end subroutine scalarPW















!subroutine my_cross(a, b, c)
!implicit none
!    real ( kind = 8 ), intent(in) :: a(3)
!	complex ( kind = 8 ), intent(in) :: b(3)
!	complex ( kind = 8 ), intent(out) :: c(3)
!        c(1) = a(2) * b(3) - a(3) * b(2)
!        c(2) = a(3) * b(1) - a(1) * b(3)
!        c(3) = a(1) * b(2) - a(2) * b(1)
!end subroutine my_cross


!subroutine my_dot(a, b, p)
!implicit none
!    real ( kind = 8 ), intent(in) :: a(3)
!	complex ( kind = 8 ), intent(in) :: b(3)
!	complex ( kind = 8 ), intent(out) :: p
!
!        p = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
!end subroutine my_dot







subroutine fieldsEDomega(omega,ep,mu,P0,points,n,E,H,vf,init)
implicit none
! a small electric dipole; that is: A=vf*exp(ima*zk*r)/r
! H=curlA
! E=1/(-ima*omega*ep)*curlcurlA
	
	!List of calling arguments
	integer *8, intent(in) :: n,init
	real ( kind = 8 ), intent(in) :: P0(3),points(12,n)
	complex ( kind = 8 ), intent(in) :: vf(3),omega,ep,mu
	complex ( kind = 8 ), intent(out) :: E(3,n),H(3,n)
	
	!List of local variables
	real ( kind = 8 ) dx,dy,dz,r
	complex ( kind = 8 ) R1,R2,au1,au2,ima,zk
	integer *8 i
	
	ima = (0.0d0,1.0d0)
	zk=omega*sqrt(ep*mu)
	do i=1,n
       if (init .eq. 0) then
			E(1,i)=0
			E(2,i)=0
			E(3,i)=0
			H(1,i)=0
			H(2,i)=0
			H(3,i)=0
		endif
			dx=points(1,i)-P0(1)
			dy=points(2,i)-P0(2)
			dz=points(3,i)-P0(3)
			r=sqrt(dx**2+dy**2+dz**2)
			R1=(ima*zk/r**2-1/r**3)
			R2=((ima*zk)**2/r**3-3*ima*zk/r**4+3/r**5)
			au1=(zk**2/r+R1)*exp(ima*zk*r)/(-ima*omega*ep)
			au2=dx*vf(1)+dy*vf(2)
			au2=au2+dz*vf(3)
			au2=au2*R2*exp(ima*zk*r)/(-ima*omega*ep)
			E(1,i)=E(1,i)+(vf(1)*au1+dx*au2)
			E(2,i)=E(2,i)+(vf(2)*au1+dy*au2)
			E(3,i)=E(3,i)+(vf(3)*au1+dz*au2)
			H(1,i)=H(1,i)+(dy*vf(3)-dz*vf(2))*R1*exp(ima*zk*r)
			H(2,i)=H(2,i)+(dz*vf(1)-dx*vf(3))*R1*exp(ima*zk*r)
			H(3,i)=H(3,i)+(dx*vf(2)-dy*vf(1))*R1*exp(ima*zk*r)
      enddo
return
end subroutine fieldsEDomega




subroutine fieldsMDomega(omega,ep,mu,P0,points,n,E,H,vf,initial)
implicit none
! a small magnetic dipole; that is: F=vf*exp(ima*zk*r)/r
! E=curlF
! H=1/(ima*omega*mu)*curlcurlF

	!List of calling arguments
	integer *8, intent(in) :: n,initial
	real ( kind = 8 ), intent(in) :: P0(3),points(12,n)
	complex ( kind = 8 ), intent(in) :: vf(3),omega,ep,mu
	complex ( kind = 8 ), intent(out) :: E(3,n),H(3,n)
	
	!List of local variables
	real ( kind = 8 ) dx,dy,dz,r,zk
	complex ( kind = 8 ) R1,R2,au1,au2,ima
	integer *8 i
	
	zk=omega*sqrt(ep*mu)
	ima = (0.0d0,1.0d0)

	do i=1,n
		if (initial .eq. 0) then
			E(1,i)=0
			E(2,i)=0
			E(3,i)=0
			H(1,i)=0
			H(2,i)=0
			H(3,i)=0
       endif
		dx=points(1,i)-P0(1)
		dy=points(2,i)-P0(2)
		dz=points(3,i)-P0(3)
		r=sqrt(dx**2+dy**2+dz**2)
		R1=(ima*zk/r**2-1/r**3)
		R2=((ima*zk)**2/r**3-3*ima*zk/r**4+3/r**5)
		au1=(zk**2/r+R1)*exp(ima*zk*r)/(ima*omega*mu)
		au2=dx*vf(1)+dy*vf(2)
		au2=au2+dz*vf(3)
		au2=au2*R2*exp(ima*zk*r)/(ima*omega*mu)
		H(1,i)=H(1,i)+(vf(1)*au1+dx*au2)
		H(2,i)=H(2,i)+(vf(2)*au1+dy*au2)
		H(3,i)=H(3,i)+(vf(3)*au1+dz*au2)
		E(1,i)=E(1,i)+(dy*vf(3)-dz*vf(2))*R1*exp(ima*zk*r)
		E(2,i)=E(2,i)+(dz*vf(1)-dx*vf(3))*R1*exp(ima*zk*r)
		E(3,i)=E(3,i)+(dx*vf(2)-dy*vf(1))*R1*exp(ima*zk*r)
	enddo
      
return
end subroutine fieldsMDomega
!
!
!
!
!

subroutine fieldsEDdpie(zk,P0,points,n,ndtarg,vf,init,E,H,A,pot,gradpot)
implicit none
! a small electric dipole; that is: A=vf*exp(ima*zk*r)/r
! H=curlA
! E=1/(-ima*zk)*curlcurlA
	
	!List of calling arguments
	integer *8, intent(in) :: n,init,ndtarg
	real ( kind = 8 ), intent(in) :: P0(3),points(ndtarg,n)
	complex ( kind = 8 ), intent(in) :: vf(3),zk
	complex ( kind = 8 ), intent(out) :: E(3,n),H(3,n),A(3,n)
	complex ( kind = 8 ), intent(out) :: pot(n),gradpot(3,n)
	
	!List of local variables
	real ( kind = 8 ) dx,dy,dz,r
	complex ( kind = 8 ) R1,R2,au1,au2,ima,myexp
	complex ( kind = 8 ) e_aux1,e_aux2,e_aux3,a_aux1,a_aux2,a_aux3
	integer *8 i
	
	ima = (0.0d0,1.0d0)

      if (init .eq. 0) then
	    do i=1,n
			E(1,i)=0.0d0
			E(2,i)=0.0d0
			E(3,i)=0.0d0
			H(1,i)=0.0d0
			H(2,i)=0.0d0
			H(3,i)=0.0d0
			A(1,i)=0.0d0
			A(2,i)=0.0d0
			A(3,i)=0.0d0
			pot(i)=0.0d0
			gradpot(1,i)=0.0d0
			gradpot(2,i)=0.0d0
			gradpot(3,i)=0.0d0
	    enddo
      endif
	  do i=1,n
    	dx=points(1,i)-P0(1)
		dy=points(2,i)-P0(2)
		dz=points(3,i)-P0(3)
		r=sqrt(dx**2+dy**2+dz**2)
		myexp=exp(ima*zk*r)
		R1=(ima*zk/r**2-1/r**3)
		R2=((ima*zk)**2/r**3-3*ima*zk/r**4+3/r**5)
		au1=(zk**2/r+R1)*myexp/(-ima*zk)
		au2=dx*vf(1)+dy*vf(2)
		au2=au2+dz*vf(3)
		au2=au2*R2*myexp/(-ima*zk)
		e_aux1=(vf(1)*au1+dx*au2)
		E(1,i)=E(1,i)+e_aux1
		e_aux2=(vf(2)*au1+dy*au2)
		E(2,i)=E(2,i)+e_aux2
		e_aux3=(vf(3)*au1+dz*au2)
		E(3,i)=E(3,i)+e_aux3
		H(1,i)=H(1,i)+(dy*vf(3)-dz*vf(2))*R1*myexp
		H(2,i)=H(2,i)+(dz*vf(1)-dx*vf(3))*R1*myexp
		H(3,i)=H(3,i)+(dx*vf(2)-dy*vf(1))*R1*myexp
		pot(i)=pot(i)+R1*myexp/(-ima*zk)*(dx*vf(1)+dy*vf(2)+dz*vf(3))
		a_aux1=myexp/r*vf(1)
		A(1,i)=A(1,i)+a_aux1
		a_aux2=myexp/r*vf(2)
		A(2,i)=A(2,i)+a_aux2
		a_aux3=myexp/r*vf(3)
		A(3,i)=A(3,i)+a_aux3
		gradpot(1,i)=gradpot(1,i)+e_aux1-ima*zk*a_aux1
		gradpot(2,i)=gradpot(2,i)+e_aux2-ima*zk*a_aux2
		gradpot(3,i)=gradpot(3,i)+e_aux3-ima*zk*a_aux3
      enddo

!      do i=1,n
!	    dx=points(1,i)-P0(1)
!		dy=points(2,i)-P0(2)
!		dz=points(3,i)-P0(3)
!		r=sqrt(dx**2+dy**2+dz**2)
!		myexp=exp(ima*zk*r)
!		R1=(ima*zk/r**2-1/r**3)

!       pot(i)=pot(i)+myexp/r
!		gradpot(1,i)=gradpot(1,i)+dx*R1*myexp
!	    gradpot(2,i)=gradpot(2,i)+dy*R1*myexp
!		gradpot(3,i)=gradpot(3,i)+dz*R1*myexp
!	  enddo


return
end subroutine fieldsEDdpie


subroutine fieldsMDdpie(zk,P0,points,n,ndtarg,vf,initial,E,H,A,pot,gradpot)
implicit none
! a small magnetic dipole; that is: F=vf*exp(ima*zk*r)/r
! E=curlF
! H=1/(ima*zk)*curlcurlF

	!List of calling arguments
	integer *8, intent(in) :: n,initial,ndtarg
	real ( kind = 8 ), intent(in) :: P0(3),points(ndtarg,n)
	complex ( kind = 8 ), intent(in) :: vf(3),zk
	complex ( kind = 8 ), intent(out) :: E(3,n),H(3,n),A(3,n)
	complex ( kind = 8 ), intent(out) :: pot(n), gradpot(3,n)
	
	!List of local variables
	real ( kind = 8 ) dx,dy,dz,r
	complex ( kind = 8 ) R1,R2,au1,au2,ima,e_aux1,e_aux2,e_aux3,myexp
	integer *8 i
	
	ima = (0.0d0,1.0d0)

    if (initial .eq. 0) then
	  do i=1,n
        E(1,i)=0.0d0
        E(2,i)=0.0d0
        E(3,i)=0.0d0
        H(1,i)=0.0d0
        H(2,i)=0.0d0
        H(3,i)=0.0d0
        A(1,i)=0.0d0
        A(2,i)=0.0d0
        A(3,i)=0.0d0
		pot(i)=0.0d0
		gradpot(1,i)=0.0d0
		gradpot(2,i)=0.0d0
		gradpot(3,i)=0.0d0
	  enddo
    endif

	do i=1,n
		dx=points(1,i)-P0(1)
		dy=points(2,i)-P0(2)
		dz=points(3,i)-P0(3)
		r=sqrt(dx**2+dy**2+dz**2)
		myexp=exp(ima*zk*r)
		R1=(ima*zk/r**2-1/r**3)
		R2=((ima*zk)**2/r**3-3*ima*zk/r**4+3/r**5)
		au1=(zk**2/r+R1)*myexp/(ima*zk)
		au2=dx*vf(1)+dy*vf(2)
		au2=au2+dz*vf(3)
		au2=au2*R2*myexp/(ima*zk)
		H(1,i)=H(1,i)+(vf(1)*au1+dx*au2)
		H(2,i)=H(2,i)+(vf(2)*au1+dy*au2)
		H(3,i)=H(3,i)+(vf(3)*au1+dz*au2)
		e_aux1=(dy*vf(3)-dz*vf(2))*R1*myexp
		E(1,i)=E(1,i)+e_aux1
		e_aux2=(dz*vf(1)-dx*vf(3))*R1*myexp
		E(2,i)=E(2,i)+e_aux2
		e_aux3=(dx*vf(2)-dy*vf(1))*R1*myexp
		E(3,i)=E(3,i)+e_aux3
		A(1,i)=A(1,i)+e_aux1/(ima*zk)
	    A(2,i)=A(2,i)+e_aux2/(ima*zk)
		A(3,i)=A(3,i)+e_aux3/(ima*zk)
!		pot(i)=pot(i)+0.0d0
!		gradpot(1,i)=gradpot(1,i)+0.0d0
!		gradpot(2,i)=gradpot(2,i)+0.0d0
!		gradpot(3,i)=gradpot(3,i)+0.0d0
	enddo
      
return
end subroutine fieldsMDdpie
