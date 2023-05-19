
subroutine kernel_Vector_Helmholtz(nk,srcinfo,targinfo,zk,E_mat)
implicit none

	!Given the representation: E=curlSk(a)+Sk(n·rho)+Sk(b)+gradSk(lambda)
	
	!And the test nx: -div: nxcurl: n· (that is the correct ordering)

    !List of calling arguments
	real ( kind = 8 ), intent(in) :: srcinfo(12)
	real ( kind = 8 ), intent(in) :: targinfo(12)
	complex ( kind = 8 ), intent(in) :: zk
    integer, intent(in) :: nk
	complex ( kind = 8 ), intent(out) :: E_mat(6,6)
	
	!List of local variables
	real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ) ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ) du(3), dv(3),sour(3),targ(3)
	real ( kind = 8 ) r, dr(3)
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),xprod_aux3(3),xprod_aux4(3)
	
	
	complex ( kind = 8 ) nxcurlSka(2,2),nxSknrho(2,1),nxSkb(2,2),nxgradSklambda(2,1)
	complex ( kind = 8 ) Dkrho(1,1),divSkb(1,2),Sklambda(1,1)
	complex ( kind = 8 ) nxcurlcurlSka(2,2),nxcurlSknrho(2,1)
	complex ( kind = 8 ) ncurlSka(1,2),nSknrho(1,1),nSkb(1,2),ngradSklambda(1,1)
	
	complex ( kind = 8 ) R1, R2, ima, my_exp

	ima=(0.0d0,1.0d0)

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
	my_exp=exp(ima*zk*r)
	
	R1=(ima*zk*r-1.0d0)/r**3*my_exp
	R2=((ima*zk)**2/r**3-3.0d0*ima*zk/r**4+3.0d0/r**5)*my_exp
	
	call orthonormalize(srcinfo(4:6),srcinfo(10:12),ru_s,rv_s)
	call orthonormalize(targinfo(4:6),targinfo(10:12),ru_t,rv_t)
	
	!nxcurlSk(a)
	call get_nxcurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,nxcurlSka)

	E_mat(1,1)=nxcurlSka(1,1)
	E_mat(1,2)=nxcurlSka(1,2)
	E_mat(2,1)=nxcurlSka(2,1)
	E_mat(2,2)=nxcurlSka(2,2)

	!nxSk(n·rho)
	call get_nxSknrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp,r,nxSknrho)
	E_mat(1,3)=nxSknrho(1,1)
	E_mat(2,3)=nxSknrho(2,1)
	
	!nxSk(b)
	call get_nxSkb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp,r,nxSkb)
	E_mat(1,4)=nxSkb(1,1)
	E_mat(1,5)=nxSkb(1,2)
	E_mat(2,4)=nxSkb(2,1)
	E_mat(2,5)=nxSkb(2,2)
	
	!nxgradSk(lambda)
	call get_nxgradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,nxgradSklambda)	
	E_mat(1,6)=nxgradSklambda(1,1)
	E_mat(2,6)=nxgradSklambda(2,1)
	
	!-div(curlSk(a))
	E_mat(3,1)=0.0d0
	E_mat(3,2)=0.0d0
	
	!-div(Sk(n·rho))=Dk(rho)
	call get_Dkrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,Dkrho)
	E_mat(3,3)=Dkrho(1,1)
	
	!-div(Sk(b))
	call get_divSkb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,divSkb)
	E_mat(3,4)=-divSkb(1,1)
	E_mat(3,5)=-divSkb(1,2)
	
	!-divgradSk(lambda)
	call get_Sklambda(my_exp,r,Sklambda)
	E_mat(3,6)=zk**2*Sklambda(1,1)
	
	!nxcurlcurlSk(a)
	call get_nxcurlcurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,R2,zk,my_exp,r,nxcurlcurlSka)
	E_mat(4,1)=nxcurlcurlSka(1,1)
	E_mat(4,2)=nxcurlcurlSka(1,2)
	E_mat(5,1)=nxcurlcurlSka(2,1)
	E_mat(5,2)=nxcurlcurlSka(2,2)
		
	!nxcurlSknrho
	call get_nxcurlSknrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,nxcurlSknrho)
	E_mat(4,3)=nxcurlSknrho(1,1)
	E_mat(5,3)=nxcurlSknrho(2,1)
	
	!nxcurlSkb
	E_mat(4,4)=nxcurlSka(1,1)
	E_mat(4,5)=nxcurlSka(1,2)	
	E_mat(5,4)=nxcurlSka(2,1)
	E_mat(5,5)=nxcurlSka(2,2)
	
	!ncurlSka
	call get_ncurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,ncurlSka)
	E_mat(6,1)=ncurlSka(1,1)
	E_mat(6,2)=ncurlSka(1,2)

	!nSknrho
	call get_nSknrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp,r,nSknrho)	
	E_mat(6,3)=nSknrho(1,1)
	
	!nSkb
	call get_nSkb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp,r,nSkb)
	E_mat(6,4)=nSkb(1,1)
	E_mat(6,5)=nSkb(1,2)
	
	!ngradSklambda
	call get_ngradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,ngradSklambda)
	E_mat(6,6)=ngradSklambda(1,1)
	
return
end subroutine kernel_Vector_Helmholtz




subroutine get_nxcurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,nxcurlSka)
implicit none

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ), intent(in) :: ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ), intent(in) :: dr(3),r
	complex ( kind = 8 ), intent(in) :: R1,my_exp
	complex ( kind = 8 ), intent(out) :: nxcurlSka(2,2)
	
 	!List of local variables
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3)

	call my_cross_v2(dr,ru_s,xprod_aux1)
	call my_cross_v2(dr,rv_s,xprod_aux2)
	
	nxcurlSka(1,1)=-DOT_PRODUCT(xprod_aux1,rv_t)*R1
	nxcurlSka(1,2)=-DOT_PRODUCT(xprod_aux2,rv_t)*R1
	nxcurlSka(2,1)=DOT_PRODUCT(xprod_aux1,ru_t)*R1
	nxcurlSka(2,2)=DOT_PRODUCT(xprod_aux2,ru_t)*R1

return
end subroutine get_nxcurlSka



subroutine get_nxSknrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp,r,nxSknrho)
implicit none
	!nxSk(n·rho)

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ), intent(in) :: ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ), intent(in) :: dr(3),r
	complex ( kind = 8 ), intent(in) :: my_exp
	complex ( kind = 8 ), intent(out) :: nxSknrho(2,1)
	
 	!List of local variables
	real ( kind = 8 ) xprod_aux1(3)

	call my_cross_v2(n_t,n_s,xprod_aux1)
	
	nxSknrho(1,1)=DOT_PRODUCT(ru_t,xprod_aux1)*my_exp/r
	nxSknrho(2,1)=DOT_PRODUCT(rv_t,xprod_aux1)*my_exp/r

return
end subroutine get_nxSknrho


subroutine get_nxSkb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp,r,nxSkb)
implicit none
	!nxSk(b)

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ), intent(in) :: ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ), intent(in) :: dr(3),r
	complex ( kind = 8 ), intent(in) :: my_exp
	complex ( kind = 8 ), intent(out) :: nxSkb(2,2)
	
 	!List of local variables
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3)

	nxSkb(1,1)=-DOT_PRODUCT(rv_t,ru_s)*my_exp/r
	nxSkb(1,2)=-DOT_PRODUCT(rv_t,rv_s)*my_exp/r
	nxSkb(2,1)=DOT_PRODUCT(ru_t,ru_s)*my_exp/r
	nxSkb(2,2)=DOT_PRODUCT(ru_t,rv_s)*my_exp/r

return
end subroutine get_nxSkb


subroutine get_nxnxSkb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp,r,nxnxSkb)
implicit none
	!nxSk(b)

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ), intent(in) :: ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ), intent(in) :: dr(3),r
	complex ( kind = 8 ), intent(in) :: my_exp
	complex ( kind = 8 ), intent(out) :: nxnxSkb(2,2)
	
 	!List of local variables
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3)

	nxnxSkb(1,1)=-DOT_PRODUCT(ru_t,ru_s)*my_exp/r
	nxnxSkb(1,2)=-DOT_PRODUCT(ru_t,rv_s)*my_exp/r
	nxnxSkb(2,1)=-DOT_PRODUCT(rv_t,ru_s)*my_exp/r
	nxnxSkb(2,2)=-DOT_PRODUCT(rv_t,rv_s)*my_exp/r

return
end subroutine get_nxnxSkb




subroutine get_nxSknxb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp,r,nxSknxb)
implicit none
	!nxSk(b)

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ), intent(in) :: ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ), intent(in) :: dr(3),r
	complex ( kind = 8 ), intent(in) :: my_exp
	complex ( kind = 8 ), intent(out) :: nxSknxb(2,2)
	
 	!List of local variables
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3)

	nxSknxb(1,1)=-DOT_PRODUCT(rv_t,rv_s)*my_exp/r
	nxSknxb(1,2)=DOT_PRODUCT(rv_t,ru_s)*my_exp/r
	nxSknxb(2,1)=DOT_PRODUCT(ru_t,rv_s)*my_exp/r
	nxSknxb(2,2)=-DOT_PRODUCT(ru_t,ru_s)*my_exp/r

return
end subroutine get_nxSknxb


subroutine get_nxgradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,nxgradSklambda)
implicit none
	!nxgradSk(lambda)

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ), intent(in) :: ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ), intent(in) :: dr(3),r
	complex ( kind = 8 ), intent(in) :: R1,my_exp
	complex ( kind = 8 ), intent(out) :: nxgradSklambda(2,1)
	
 	!List of local variables
	real ( kind = 8 ) xprod_aux1(3)

	!nxgradSk(lambda)
	nxgradSklambda(1,1)=-DOT_PRODUCT(rv_t,dr)*R1
	nxgradSklambda(2,1)=DOT_PRODUCT(ru_t,dr)*R1


return
end subroutine get_nxgradSklambda



subroutine get_nxnxgradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,nxnxgradSklambda)
implicit none
	!nxgradSk(lambda)

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ), intent(in) :: ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ), intent(in) :: dr(3),r
	complex ( kind = 8 ), intent(in) :: R1,my_exp
	complex ( kind = 8 ), intent(out) :: nxnxgradSklambda(2,1)
	
 	!List of local variables
	real ( kind = 8 ) xprod_aux1(3)

	!nxgradSk(lambda)
	nxnxgradSklambda(1,1)=-DOT_PRODUCT(ru_t,dr)*R1
	nxnxgradSklambda(2,1)=-DOT_PRODUCT(rv_t,dr)*R1


return
end subroutine get_nxnxgradSklambda


subroutine get_Dkrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,Dkrho)
implicit none
	!-div(Sk(n·rho))=Dk(rho)

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ), intent(in) :: ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ), intent(in) :: dr(3),r
	complex ( kind = 8 ), intent(in) :: R1,my_exp
	complex ( kind = 8 ), intent(out) :: Dkrho(1,1)
	

	!-div(Sk(n·rho))
	Dkrho(1,1)=-DOT_PRODUCT(n_s,dr)*R1
	
return
end subroutine get_Dkrho


subroutine get_divSkb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,divSkb)
implicit none
	!divSk(b)

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ), intent(in) :: ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ), intent(in) :: dr(3),r
	complex ( kind = 8 ), intent(in) :: R1,my_exp
	complex ( kind = 8 ), intent(out) :: divSkb(1,2)
	
 	!List of local variables
	real ( kind = 8 ) xprod_aux1(3)

	!div(Sk(b))
	divSkb(1,1)=DOT_PRODUCT(dr,ru_s)*R1
	divSkb(1,2)=DOT_PRODUCT(dr,rv_s)*R1

return
end subroutine get_divSkb



subroutine get_divSknxb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,divSknxb)
implicit none
	!divSk(b)

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ), intent(in) :: ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ), intent(in) :: dr(3),r
	complex ( kind = 8 ), intent(in) :: R1,my_exp
	complex ( kind = 8 ), intent(out) :: divSknxb(1,2)
	
 	!List of local variables
	real ( kind = 8 ) xprod_aux1(3)

	!div(Sk(b))
	divSknxb(1,1)=DOT_PRODUCT(dr,rv_s)*R1
	divSknxb(1,2)=-DOT_PRODUCT(dr,ru_s)*R1

return
end subroutine get_divSknxb


subroutine get_Sklambda(my_exp,r,Sklambda)
implicit none
	!Sk(lambda)

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: r
	complex ( kind = 8 ), intent(in) :: my_exp
	complex ( kind = 8 ), intent(out) :: Sklambda(1,1)
	
	Sklambda(1,1)=my_exp/r
	
return
end subroutine get_Sklambda


subroutine get_nxcurlcurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,R2,zk,my_exp,r,nxcurlcurlSka)
implicit none
	!nxcurlcurlSk(a)

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ), intent(in) :: ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ), intent(in) :: dr(3),r
	complex ( kind = 8 ), intent(in) :: R1,R2,zk,my_exp
	complex ( kind = 8 ), intent(out) :: nxcurlcurlSka(2,2)
	
 	!List of local variables
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3)
	complex ( kind = 8 ) nxSkb(2,2)	
	
	nxSkb(1,1)=-DOT_PRODUCT(rv_t,ru_s)*my_exp/r
	nxSkb(1,2)=-DOT_PRODUCT(rv_t,rv_s)*my_exp/r
	nxSkb(2,1)=DOT_PRODUCT(ru_t,ru_s)*my_exp/r
	nxSkb(2,2)=DOT_PRODUCT(ru_t,rv_s)*my_exp/r

	nxcurlcurlSka(1,1)=zk**2*nxSkb(1,1)
	nxcurlcurlSka(1,2)=zk**2*nxSkb(1,2)
	nxcurlcurlSka(2,1)=zk**2*nxSkb(2,1)
	nxcurlcurlSka(2,2)=zk**2*nxSkb(2,2)
	
	nxcurlcurlSka(1,1)=nxcurlcurlSka(1,1)+(-DOT_PRODUCT(rv_t,ru_s)*R1)
	nxcurlcurlSka(1,2)=nxcurlcurlSka(1,2)+(-DOT_PRODUCT(rv_t,rv_s)*R1)
	nxcurlcurlSka(2,1)=nxcurlcurlSka(2,1)+(+DOT_PRODUCT(ru_t,ru_s)*R1)
	nxcurlcurlSka(2,2)=nxcurlcurlSka(2,2)+(+DOT_PRODUCT(ru_t,rv_s)*R1)
	
	nxcurlcurlSka(1,1)=nxcurlcurlSka(1,1)-(-DOT_PRODUCT(rv_t,dr)*DOT_PRODUCT(ru_s,-dr)*R2)
	nxcurlcurlSka(1,2)=nxcurlcurlSka(1,2)-(-DOT_PRODUCT(rv_t,dr)*DOT_PRODUCT(rv_s,-dr)*R2)
	nxcurlcurlSka(2,1)=nxcurlcurlSka(2,1)-(+DOT_PRODUCT(ru_t,dr)*DOT_PRODUCT(ru_s,-dr)*R2)
	nxcurlcurlSka(2,2)=nxcurlcurlSka(2,2)-(+DOT_PRODUCT(ru_t,dr)*DOT_PRODUCT(rv_s,-dr)*R2)

return
end subroutine get_nxcurlcurlSka


subroutine get_ngradDkrho(n_s,n_t,dr,R1,R2,zk,my_exp,r,ngradDkrho)
implicit none
	!nxcurlcurlSk(a)

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: n_s(3)
	real ( kind = 8 ), intent(in) :: n_t(3)
	real ( kind = 8 ), intent(in) :: dr(3),r
	complex ( kind = 8 ), intent(in) :: R1,R2,zk,my_exp
	complex ( kind = 8 ), intent(out) :: ngradDkrho(1,1)
	
 	!List of local variables
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3)
	complex ( kind = 8 ) nxSkb(2,2)	
	
	ngradDkrho(1,1)=DOT_PRODUCT(n_t,dr)*DOT_PRODUCT(n_s,-dr)*R2
	ngradDkrho(1,1)=ngradDkrho(1,1)-DOT_PRODUCT(n_t,n_s)*R1
	
return
end subroutine get_ngradDkrho

subroutine get_nxcurlSknrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,nxcurlSknrho)
implicit none
	!nxgradSk(lambda)

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ), intent(in) :: ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ), intent(in) :: dr(3),r
	complex ( kind = 8 ), intent(in) :: R1,my_exp
	complex ( kind = 8 ), intent(out) :: nxcurlSknrho(2,1)
	
 	!List of local variables
	real ( kind = 8 ) xprod_aux1(3)

	call my_cross_v2(dr, n_s, xprod_aux1)
	nxcurlSknrho(1,1)=-DOT_PRODUCT(rv_t,xprod_aux1)*R1
	nxcurlSknrho(2,1)=DOT_PRODUCT(ru_t,xprod_aux1)*R1

return
end subroutine get_nxcurlSknrho


subroutine get_ncurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,ncurlSka)
implicit none
	!ncurnSk(a)

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ), intent(in) :: ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ), intent(in) :: dr(3),r
	complex ( kind = 8 ), intent(in) :: R1,my_exp
	complex ( kind = 8 ), intent(out) :: ncurlSka(1,2)
	
 	!List of local variables
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3)

	call my_cross_v2(dr, ru_s, xprod_aux1)
	call my_cross_v2(dr, rv_s, xprod_aux2)
	ncurlSka(1,1)=DOT_PRODUCT(n_t,xprod_aux1)*R1
	ncurlSka(1,2)=DOT_PRODUCT(n_t,xprod_aux2)*R1

return
end subroutine get_ncurlSka


subroutine get_nSknrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp,r,nSknrho)
implicit none
	!nSk(nrho)

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ), intent(in) :: ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ), intent(in) :: dr(3),r
	complex ( kind = 8 ), intent(in) :: my_exp
	complex ( kind = 8 ), intent(out) :: nSknrho(1,1)
	
 	!List of local variables
	real ( kind = 8 ) xprod_aux1(3)

	nSknrho(1,1)=DOT_PRODUCT(n_t,n_s)*my_exp/r

return
end subroutine get_nSknrho



subroutine get_nSkb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp,r,nSkb)
implicit none
	!nSk(b)

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ), intent(in) :: ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ), intent(in) :: dr(3),r
	complex ( kind = 8 ), intent(in) :: my_exp
	complex ( kind = 8 ), intent(out) :: nSkb(1,2)
	
 	!List of local variables
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3)

	nSkb(1,1)=DOT_PRODUCT(n_t,ru_s)*my_exp/r
	nSkb(1,2)=DOT_PRODUCT(n_t,rv_s)*my_exp/r

return
end subroutine get_nSkb
!
!
!
!
!
subroutine get_nSknxb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp,r,nSknxb)
implicit none
	!nSk(nxb)

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ), intent(in) :: ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ), intent(in) :: dr(3),r
	complex ( kind = 8 ), intent(in) :: my_exp
	complex ( kind = 8 ), intent(out) :: nSknxb(1,2)
	
 	!List of local variables
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3)

	nSknxb(1,1)=DOT_PRODUCT(n_t,rv_s)*my_exp/r
	nSknxb(1,2)=-DOT_PRODUCT(n_t,ru_s)*my_exp/r

return
end subroutine get_nSknxb


subroutine get_ngradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,ngradSklambda)
implicit none
	!ngradSk(lambda)

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ), intent(in) :: ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ), intent(in) :: dr(3),r
	complex ( kind = 8 ), intent(in) :: R1,my_exp
	complex ( kind = 8 ), intent(out) :: ngradSklambda(1,1)
	
	ngradSklambda(1,1)=DOT_PRODUCT(n_t,dr)*R1
	
return
end subroutine get_ngradSklambda


! subroutine orthonormalize_all(du,normal,ru,rv,ns)
! !f2py intent(in) du, normal, ns
! !f2py intent(out) ru,rv
! implicit none

!     !List of calling arguments
! 	integer, intent(in) :: ns
! 	real ( kind = 8 ), intent(in) :: du(3,ns), normal(3,ns)
! 	real ( kind = 8 ), intent(out) :: ru(3,ns), rv(3,ns)

! 	!List of local variables
! 	real ( kind = 8 ) aux
! 	integer count1

! 	do count1=1,ns
! 		call orthonormalize(du(:,count1),normal(:,count1),ru(:,count1),rv(:,count1))
! 	enddo

! return
! end subroutine orthonormalize_all






! subroutine orthonormalize(du,normal,ru,rv)
!   implicit none

!   !List of calling arguments
!   real ( kind = 8 ), intent(in) :: du(3), normal(3)
!   real ( kind = 8 ), intent(out) :: ru(3), rv(3)

!   !List of local variables
!   real ( kind = 8 ) aux

!   aux=sqrt(du(1)**2+du(2)**2+du(3)**2)

!   ru(1)=du(1)/aux
!   ru(2)=du(2)/aux
!   ru(3)=du(3)/aux

!   call my_cross_v2(normal, ru, rv)
!   return
! end subroutine orthonormalize






 subroutine my_cross_v2(a, b, c)
 implicit none

     real ( kind = 8 ), intent(in) :: a(3),b(3)
     real ( kind = 8 ), intent(out) :: c(3)

         c(1) = a(2) * b(3) - a(3) * b(2)
         c(2) = a(3) * b(1) - a(1) * b(3)
         c(3) = a(1) * b(2) - a(2) * b(1)

 end subroutine my_cross_v2
