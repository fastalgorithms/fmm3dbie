



subroutine rwg2vr(srcvals,ns,norder,ntri,npoints,points,triang,rwg,nrwg,coefrwg,coefvr)
implicit none

  !List of calling arguments
  integer, intent(in) :: ntri,npoints,ns,nrwg,norder
  integer, intent(in) :: triang(3,ntri)
  real ( kind = 8 ), intent(in) :: points(3,npoints),srcvals(12,ns)
  integer, intent(in) :: rwg(4,nrwg)
  complex ( kind = 8 ), intent(in) :: coefrwg(nrwg)
  complex ( kind = 8 ), intent(out) :: coefvr(3*ns)

  !List of local variables
  integer count1,count2,count3,count4,flag,nspace,icount,sol,t1,t2,caux
  integer npols,i1,i2,v1i,v2i,v3i
  real ( kind = 8 ) P1(3),P2(3),P3(3),r_aux,A_t,ln,n_vect(3)
  real ( kind = 8 ) Pa(3),Pb(3),Pc(3)
  real ( kind = 8 ), allocatable :: rx(:),ry(:),rz(:)
	real ( kind = 8 ), allocatable :: uv(:,:),umatr(:,:),vmatr(:,:),w(:)
  real ( kind = 8 ), allocatable :: u_vect_s(:,:),v_vect_s(:,:)

  do count1=1,3*ns
    coefvr(count1)=0.0d0
  enddo

  allocate(u_vect_s(3,ns),v_vect_s(3,ns))

  call orthonormalize_all_v2(srcvals(4:6,:),srcvals(10:12,:),&
  &u_vect_s,v_vect_s,ns)
  
  npols = (norder+1)*(norder+2)/2
  allocate(rx(npols),ry(npols),rz(npols))

  allocate(uv(2,npols),umatr(npols,npols),vmatr(npols,npols),w(npols))

  call vioreanu_simplex_quad(norder,npols,uv,umatr,vmatr,w)

  do count1=1,nrwg
    t1=rwg(1,count1)
    v1i=rwg(3,count1)
    v2i=v1i+1
    v3i=v1i+2
    if (v2i>3) then
      v2i=v2i-3
    endif
    if (v3i>3) then
      v3i=v3i-3
    endif
    P1=points(:,triang(v1i,t1))
    P2=points(:,triang(v2i,t1))
    P3=points(:,triang(v3i,t1))

!    write (*,*) 'P1: ',P1

    Pa(:)=points(:,triang(1,t1))
    Pb(:)=points(:,triang(2,t1))
    Pc(:)=points(:,triang(3,t1))

    rx=(Pb(1)-Pa(1))*uv(1,:)+(Pc(1)-Pa(1))*uv(2,:)+Pa(1)-P1(1)
    ry=(Pb(2)-Pa(2))*uv(1,:)+(Pc(2)-Pa(2))*uv(2,:)+Pa(2)-P1(2)
    rz=(Pb(3)-Pa(3))*uv(1,:)+(Pc(3)-Pa(3))*uv(2,:)+Pa(3)-P1(3)

    ! ahora falta escalarlo y proyectar
    call my_cross_v3(P2-P1, P3-P1, n_vect)
    A_t=norm2(n_vect)/2.0d0
    ln=norm2(P3-P2)
    r_aux=ln/(2.0d0*A_t)
    rx=rx*r_aux
    ry=ry*r_aux
    rz=rz*r_aux
    i1=npols*(t1-1)+1
    i2=npols*t1
    coefvr(i1:i2)=coefvr(i1:i2)+(rx*u_vect_s(1,i1:i2)+&
     &ry*u_vect_s(2,i1:i2)+rz*u_vect_s(3,i1:i2))*coefrwg(count1)
    coefvr((ns+i1):(ns+i2))=coefvr((ns+i1):(ns+i2))+(rx*v_vect_s(1,i1:i2)&
     &+ry*v_vect_s(2,i1:i2)+rz*v_vect_s(3,i1:i2))*coefrwg(count1)

    coefvr((2*ns+i1):(2*ns+i2))=coefvr((2*ns+i1):(2*ns+i2))+2.0d0*r_aux*coefrwg(count1)


    t1=rwg(2,count1)
    v1i=rwg(4,count1)
    v2i=v1i+1
    v3i=v1i+2
    if (v2i>3) then
      v2i=v2i-3
    endif
    if (v3i>3) then
      v3i=v3i-3
    endif
    P1=points(:,triang(v1i,t1))
    P2=points(:,triang(v2i,t1))
    P3=points(:,triang(v3i,t1))

!    write (*,*) 'P1: ',P1

    Pa(:)=points(:,triang(1,t1))
    Pb(:)=points(:,triang(2,t1))
    Pc(:)=points(:,triang(3,t1))

    rx=(Pb(1)-Pa(1))*uv(1,:)+(Pc(1)-Pa(1))*uv(2,:)+Pa(1)-P1(1)
    ry=(Pb(2)-Pa(2))*uv(1,:)+(Pc(2)-Pa(2))*uv(2,:)+Pa(2)-P1(2)
    rz=(Pb(3)-Pa(3))*uv(1,:)+(Pc(3)-Pa(3))*uv(2,:)+Pa(3)-P1(3)

    ! ahora falta escalarlo y proyectar
    call my_cross_v3(P2-P1, P3-P1, n_vect)
    A_t=norm2(n_vect)/2.0d0
    ln=norm2(P3-P2)
    r_aux=ln/(2.0d0*A_t)
    rx=-rx*r_aux
    ry=-ry*r_aux
    rz=-rz*r_aux
    i1=npols*(t1-1)+1
    i2=npols*t1
    coefvr(i1:i2)=coefvr(i1:i2)+(rx*u_vect_s(1,i1:i2)+&
     &ry*u_vect_s(2,i1:i2)+rz*u_vect_s(3,i1:i2))*coefrwg(count1)
    coefvr((ns+i1):(ns+i2))=coefvr((ns+i1):(ns+i2))+(rx*v_vect_s(1,i1:i2)+&
     &ry*v_vect_s(2,i1:i2)+rz*v_vect_s(3,i1:i2))*coefrwg(count1)
    coefvr((2*ns+i1):(2*ns+i2))=coefvr((2*ns+i1):(2*ns+i2))-2.0d0*r_aux*coefrwg(count1)

  enddo
!  write (*,*) 'debug rwg: '
!  do count1=1,ns
!    write (*,*) u_vect_s(:,count1),v_vect_s(:,count1),real(coefvr(count1)),real(coefvr(count1+ns)),real(coefvr(count1+2*ns))
!  enddo
!  read (*,*)
!  stop

return
end subroutine rwg2vr






subroutine vr2rwg_vector(srcvals,wts,ns,norder,ntri,npoints,&
    points,triang,rwg,nrwg,coefvr,coefrwg)
implicit none

  !List of calling arguments
  integer, intent(in) :: ntri,npoints,ns,nrwg,norder
  integer, intent(in) :: triang(3,ntri)
  real ( kind = 8 ), intent(in) :: points(3,npoints),srcvals(12,ns),wts(ns)
  integer, intent(in) :: rwg(4,nrwg)
  complex ( kind = 8 ), intent(out) :: coefrwg(nrwg)
  complex ( kind = 8 ), intent(in) :: coefvr(2*ns)

  !List of local variables
  integer count1,count2,count3,count4,flag,nspace,icount,sol,t1,t2,caux
  integer npols,i1,i2,v1i,v2i,v3i
  real ( kind = 8 ) P1(3),P2(3),P3(3),r_aux,A_t,ln,n_vect(3)
  real ( kind = 8 ) Pa(3),Pb(3),Pc(3)
  real ( kind = 8 ), allocatable :: rx(:),ry(:),rz(:)
	real ( kind = 8 ), allocatable :: uv(:,:),umatr(:,:),vmatr(:,:),w(:)
  real ( kind = 8 ), allocatable :: u_vect_s(:,:),v_vect_s(:,:)
  complex ( kind = 8 ), allocatable :: E(:,:)

  do count1=1,nrwg
    coefrwg(count1)=0.0d0
  enddo

  allocate(u_vect_s(3,ns),v_vect_s(3,ns))
  allocate(E(3,ns))

  call orthonormalize_all_v2(srcvals(4:6,:),srcvals(10:12,:),&
  &u_vect_s,v_vect_s,ns)

  npols = (norder+1)*(norder+2)/2
  allocate(rx(npols),ry(npols),rz(npols))

	allocate(uv(2,npols),umatr(npols,npols),vmatr(npols,npols),w(npols))

	call vioreanu_simplex_quad(norder,npols,uv,umatr,vmatr,w)

  do count1=1,nrwg
    t1=rwg(1,count1)
    v1i=rwg(3,count1)
    v2i=v1i+1
    v3i=v1i+2
    if (v2i>3) then
      v2i=v2i-3
    endif
    if (v3i>3) then
      v3i=v3i-3
    endif
    P1=points(:,triang(v1i,t1))
    P2=points(:,triang(v2i,t1))
    P3=points(:,triang(v3i,t1))

!    write (*,*) 'P1: ',P1

    Pa(:)=points(:,triang(1,t1))
    Pb(:)=points(:,triang(2,t1))
    Pc(:)=points(:,triang(3,t1))

    rx=(Pb(1)-Pa(1))*uv(1,:)+(Pc(1)-Pa(1))*uv(2,:)+Pa(1)-P1(1)
    ry=(Pb(2)-Pa(2))*uv(1,:)+(Pc(2)-Pa(2))*uv(2,:)+Pa(2)-P1(2)
    rz=(Pb(3)-Pa(3))*uv(1,:)+(Pc(3)-Pa(3))*uv(2,:)+Pa(3)-P1(3)

    ! ahora falta escalarlo y proyectar
    call my_cross_v3(P2-P1, P3-P1, n_vect)
    A_t=norm2(n_vect)/2.0d0
    ln=norm2(P3-P2)
    r_aux=ln/(2.0d0*A_t)
    rx=rx*r_aux
    ry=ry*r_aux
    rz=rz*r_aux
    i1=npols*(t1-1)+1
    i2=npols*t1

    E(1,i1:i2)=coefvr(i1:i2)*u_vect_s(1,i1:i2)+coefvr(ns+i1:ns+i2)*v_vect_s(1,i1:i2)
    E(2,i1:i2)=coefvr(i1:i2)*u_vect_s(2,i1:i2)+coefvr(ns+i1:ns+i2)*v_vect_s(2,i1:i2)
    E(3,i1:i2)=coefvr(i1:i2)*u_vect_s(3,i1:i2)+coefvr(ns+i1:ns+i2)*v_vect_s(3,i1:i2)

    coefrwg(count1)=coefrwg(count1)+sum((E(1,i1:i2)*rx+E(2,i1:i2)*ry+E(3,i1:i2)*rz)*wts(i1:i2))

    t1=rwg(2,count1)
    v1i=rwg(4,count1)
    v2i=v1i+1
    v3i=v1i+2
    if (v2i>3) then
      v2i=v2i-3
    endif
    if (v3i>3) then
      v3i=v3i-3
    endif
    P1=points(:,triang(v1i,t1))
    P2=points(:,triang(v2i,t1))
    P3=points(:,triang(v3i,t1))

!    write (*,*) 'P1: ',P1

    Pa(:)=points(:,triang(1,t1))
    Pb(:)=points(:,triang(2,t1))
    Pc(:)=points(:,triang(3,t1))

    rx=(Pb(1)-Pa(1))*uv(1,:)+(Pc(1)-Pa(1))*uv(2,:)+Pa(1)-P1(1)
    ry=(Pb(2)-Pa(2))*uv(1,:)+(Pc(2)-Pa(2))*uv(2,:)+Pa(2)-P1(2)
    rz=(Pb(3)-Pa(3))*uv(1,:)+(Pc(3)-Pa(3))*uv(2,:)+Pa(3)-P1(3)

    ! ahora falta escalarlo y proyectar
    call my_cross_v3(P2-P1, P3-P1, n_vect)
    A_t=norm2(n_vect)/2.0d0
    ln=norm2(P3-P2)
    r_aux=ln/(2.0d0*A_t)
    rx=-rx*r_aux
    ry=-ry*r_aux
    rz=-rz*r_aux
    i1=npols*(t1-1)+1
    i2=npols*t1

    E(1,i1:i2)=coefvr(i1:i2)*u_vect_s(1,i1:i2)+coefvr(ns+i1:ns+i2)*v_vect_s(1,i1:i2)
    E(2,i1:i2)=coefvr(i1:i2)*u_vect_s(2,i1:i2)+coefvr(ns+i1:ns+i2)*v_vect_s(2,i1:i2)
    E(3,i1:i2)=coefvr(i1:i2)*u_vect_s(3,i1:i2)+coefvr(ns+i1:ns+i2)*v_vect_s(3,i1:i2)

    coefrwg(count1)=coefrwg(count1)+sum((E(1,i1:i2)*rx+E(2,i1:i2)*ry+E(3,i1:i2)*rz)*wts(i1:i2))

  enddo

return
end subroutine vr2rwg_vector






subroutine vr2rwg_scalar(srcvals,wts,ns,norder,ntri,npoints,points,triang, &
      rwg,nrwg,coefvr,coefrwg)
implicit none

  !List of calling arguments
  integer, intent(in) :: ntri,npoints,ns,nrwg,norder
  integer, intent(in) :: triang(3,ntri)
  real ( kind = 8 ), intent(in) :: points(3,npoints),srcvals(12,ns),wts(ns)
  integer, intent(in) :: rwg(4,nrwg)
  complex ( kind = 8 ), intent(out) :: coefrwg(nrwg)
  complex ( kind = 8 ), intent(in) :: coefvr(ns)

  !List of local variables
  integer count1,count2,count3,count4,flag,nspace,icount,sol,t1,t2,caux
  integer npols,i1,i2,v1i,v2i,v3i
  real ( kind = 8 ) P1(3),P2(3),P3(3),r_aux,A_t,ln,n_vect(3)
  real ( kind = 8 ) Pa(3),Pb(3),Pc(3)
  real ( kind = 8 ), allocatable :: rx(:),ry(:),rz(:)
	real ( kind = 8 ), allocatable :: uv(:,:),umatr(:,:),vmatr(:,:),w(:)
  real ( kind = 8 ), allocatable :: u_vect_s(:,:),v_vect_s(:,:)
  complex ( kind = 8 ), allocatable :: E(:,:)

  do count1=1,nrwg
    coefrwg(count1)=0.0d0
  enddo

  allocate(u_vect_s(3,ns),v_vect_s(3,ns))
  allocate(E(3,ns))

  call orthonormalize_all_v2(srcvals(4:6,:),srcvals(10:12,:),&
  &u_vect_s,v_vect_s,ns)
  
  npols = (norder+1)*(norder+2)/2
  allocate(rx(npols),ry(npols),rz(npols))

	allocate(uv(2,npols),umatr(npols,npols),vmatr(npols,npols),w(npols))

	call vioreanu_simplex_quad(norder,npols,uv,umatr,vmatr,w)

  do count1=1,nrwg
    t1=rwg(1,count1)
    v1i=rwg(3,count1)
    v2i=v1i+1
    v3i=v1i+2
    if (v2i>3) then
      v2i=v2i-3
    endif
    if (v3i>3) then
      v3i=v3i-3
    endif
    P1=points(:,triang(v1i,t1))
    P2=points(:,triang(v2i,t1))
    P3=points(:,triang(v3i,t1))

!    write (*,*) 'P1: ',P1

    Pa(:)=points(:,triang(1,t1))
    Pb(:)=points(:,triang(2,t1))
    Pc(:)=points(:,triang(3,t1))

    rx=(Pb(1)-Pa(1))*uv(1,:)+(Pc(1)-Pa(1))*uv(2,:)+Pa(1)-P1(1)
    ry=(Pb(2)-Pa(2))*uv(1,:)+(Pc(2)-Pa(2))*uv(2,:)+Pa(2)-P1(2)
    rz=(Pb(3)-Pa(3))*uv(1,:)+(Pc(3)-Pa(3))*uv(2,:)+Pa(3)-P1(3)

    ! ahora falta escalarlo y proyectar
    call my_cross_v3(P2-P1, P3-P1, n_vect)
    A_t=norm2(n_vect)/2.0d0
    ln=norm2(P3-P2)
    r_aux=ln/(2.0d0*A_t)
    rx=rx*r_aux
    ry=ry*r_aux
    rz=rz*r_aux
    i1=npols*(t1-1)+1
    i2=npols*t1

    coefrwg(count1)=coefrwg(count1)+2.0d0*r_aux*sum(coefvr(i1:i2)*wts(i1:i2))

    t1=rwg(2,count1)
    v1i=rwg(4,count1)
    v2i=v1i+1
    v3i=v1i+2
    if (v2i>3) then
      v2i=v2i-3
    endif
    if (v3i>3) then
      v3i=v3i-3
    endif
    P1=points(:,triang(v1i,t1))
    P2=points(:,triang(v2i,t1))
    P3=points(:,triang(v3i,t1))

!    write (*,*) 'P1: ',P1

    Pa(:)=points(:,triang(1,t1))
    Pb(:)=points(:,triang(2,t1))
    Pc(:)=points(:,triang(3,t1))

    rx=(Pb(1)-Pa(1))*uv(1,:)+(Pc(1)-Pa(1))*uv(2,:)+Pa(1)-P1(1)
    ry=(Pb(2)-Pa(2))*uv(1,:)+(Pc(2)-Pa(2))*uv(2,:)+Pa(2)-P1(2)
    rz=(Pb(3)-Pa(3))*uv(1,:)+(Pc(3)-Pa(3))*uv(2,:)+Pa(3)-P1(3)

    ! ahora falta escalarlo y proyectar
    call my_cross_v3(P2-P1, P3-P1, n_vect)
    A_t=norm2(n_vect)/2.0d0
    ln=norm2(P3-P2)
    r_aux=ln/(2.0d0*A_t)
    rx=-rx*r_aux
    ry=-ry*r_aux
    rz=-rz*r_aux
    i1=npols*(t1-1)+1
    i2=npols*t1

    coefrwg(count1)=coefrwg(count1)-2.0d0*r_aux*sum(coefvr(i1:i2)*wts(i1:i2))

  enddo

return
end subroutine vr2rwg_scalar







subroutine rwg_basis(ntri,npoints,points,triang,rwg,nrwg)
implicit none

  !List of calling arguments
  integer, intent(in) :: ntri,npoints
  integer, intent(in) :: triang(3,ntri)
  real ( kind = 8 ), intent(in) :: points(3,npoints)
  integer, intent(out) :: nrwg,rwg(4,3*ntri)

  !List of local variables
  integer count1,count2,count3,count4,flag,nspace,icount,sol,t1,t2,caux
  integer, allocatable :: tri_map(:,:)
  integer, allocatable :: rwg_aux(:,:),flags(:,:)
  real ( kind = 8 ) eps

    allocate(rwg_aux(3,ntri),flags(3,ntri))

    sol=0
    nspace=20
    eps=1.0d-15
    allocate(tri_map(nspace,npoints))

    do count1=1,ntri
      rwg_aux(1,count1)=0
      rwg_aux(2,count1)=0
      rwg_aux(3,count1)=0
    enddo
    do count1=1,3*ntri
      rwg(:,count1)=0
    enddo



    do count1=1,npoints
      do count2=1,nspace
        tri_map(count2,count1)=0
      enddo
    enddo
    
    do count1=1,ntri
      call add_tri(count1,triang(1,count1),tri_map,nspace,npoints)
      call add_tri(count1,triang(2,count1),tri_map,nspace,npoints)
      call add_tri(count1,triang(3,count1),tri_map,nspace,npoints)
    enddo

!    do count1=1,npoints
!      write (*,*) count1,'tri_map',tri_map(:,count1)
!    enddo

    do count1=1,npoints
      do count2=1,nspace
        if (tri_map(count2,count1).eq.0) then
          exit
        else
          do count3=1,nspace
            if (tri_map(count3,count1).eq.0) then
              exit
            else
              t1=tri_map(count2,count1)
              t2=tri_map(count3,count1)
              call share_wedge(t1,t2,triang,ntri,sol)
              if (sol.eq.1) then
                call add_rwg(t1,t2,rwg_aux,ntri)
              endif
            endif
          enddo
        endif
      enddo
    enddo

    do count1=1,ntri
      flags(1,count1)=0
      flags(2,count1)=0
      flags(3,count1)=0
    enddo


    do count1=1,ntri
      if (rwg_aux(1,count1).ne.0) then
        flags(1,count1)=1
      endif
      if (rwg_aux(2,count1).ne.0) then
        flags(2,count1)=1
      endif
      if (rwg_aux(3,count1).ne.0) then
        flags(3,count1)=1
      endif
    enddo

    icount=1
    do count1=1,ntri
     do count2=1,3
!       write (*,*) ((rwg_aux(count2,count1).ne.0).and.(flags(count2,count1).ne.0)),ntri
       if ((rwg_aux(count2,count1).ne.0).and.(flags(count2,count1).ne.0)) then
         if (rwg_aux(1,rwg_aux(count2,count1)).eq.count1) then
           caux=1
         elseif (rwg_aux(2,rwg_aux(count2,count1)).eq.count1) then
           caux=2
         elseif (rwg_aux(3,rwg_aux(count2,count1)).eq.count1) then
           caux=3
         else
          write (*,*) 'error!!!'
          stop
         endif
         t1=count1
         t2=rwg_aux(count2,count1)
         rwg(1,icount)=t1
         rwg(2,icount)=t2
         flags(count2,count1)=0
         flags(caux,rwg_aux(count2,count1))=0
         do count3=1,3
           if ((triang(count3,t1).ne.triang(1,t2)).and.(triang(count3,t1).ne.&
            &triang(2,t2)).and.(triang(count3,t1).ne.triang(3,t2))) then
             rwg(3,icount)=count3
           endif
           if ((triang(count3,t2).ne.triang(1,t1)).and.(triang(count3,t2).ne.&
            &triang(2,t1)).and.(triang(count3,t2).ne.triang(3,t1))) then
             rwg(4,icount)=count3
           endif
         enddo
         icount=icount+1
!         do count4=1,20
!           write (*,*) count4,rwg_aux(:,count4),flags(:,count4)
!           write (*,*) count4,rwg(:,count1)
!         enddo
!         read (*,*)
       endif
     enddo
    enddo
    nrwg=icount-1
!    write (*,*) 'datoo',nrwg,3*ntri
!    do count1=1,nrwg+10
!      write (*,*) count1,rwg(:,count1)
!    enddo
return
end subroutine rwg_basis







subroutine add_rwg(t1,t2,rwg,ntri)
implicit none

  !List of calling arguments
  integer, intent(in) :: ntri,t1,t2
  integer, intent(inout) :: rwg(3,ntri)
  
  !List of local variables
  integer count1,count2,icount

  do count1=1,3
    if (rwg(count1,t1).eq.t2) then
      exit
    elseif (rwg(count1,t1).eq.0) then
      rwg(count1,t1)=t2
      exit
    endif
  enddo

  do count1=1,3
    if (rwg(count1,t2).eq.t1) then
      exit
    elseif (rwg(count1,t2).eq.0) then
      rwg(count1,t2)=t1
      exit
    endif
  enddo

!  do count1=1,ntri
!    write (*,*) rwg(:,count1)
!  enddo
!  write (*,*) t1,t2
!  read (*,*)

return
end subroutine add_rwg






subroutine share_wedge(t1,t2,triang,ntri,sol)
implicit none

  !List of calling arguments
  integer, intent(in) :: ntri,t1,t2
  integer, intent(in) :: triang(3,ntri)
  integer, intent(out) :: sol
  
  !List of local variables
  integer count1,count2,icount

    icount=0
    do count1=1,3
      do count2=1,3
        if ((triang(count1,t1)-triang(count2,t2)).eq.0) then 
          icount=icount+1
        endif
      enddo
    enddo

    if (icount.eq.2) then
      sol=1
    else
      sol=0
    endif
!    write (*,*) triang(:,t1)
!    write (*,*) triang(:,t2)
!    write (*,*) sol
!    read (*,*)

return
end subroutine share_wedge








subroutine fix_tri_geom(ntri,npoints,points,triang,npoints_out,points_out,triang_out)
implicit none

  !List of calling arguments
  integer, intent(in) :: ntri,npoints
  integer, intent(out) :: triang(3,ntri),npoints_out
  real ( kind = 8 ), intent(in) :: points(3,npoints)
  real ( kind = 8 ), intent(out) :: points_out(3,npoints)
  integer, intent(out) :: triang_out(3,ntri)


  !List of local variables
  integer count1,count2,count3,flag,nspace,icount,n_aux
  integer, allocatable :: flags(:),tri_map(:,:)
  real ( kind = 8 ), allocatable :: diam_tri(:)
  real ( kind = 8 ) eps,Pa(3),Pb(3),Pc(3),d_vect(3),d,tol
  real ( kind = 8 ) Pa1(3),Pb1(3),Pc1(3),Pa2(3),Pb2(3),Pc2(3),n1(3),n2(3)

    eps=1.0d-1
    eps=eps**2
    nspace=20   ! at most one point belongs to nspace triangles..
    allocate(flags(npoints),tri_map(nspace,npoints))
    allocate(diam_tri(ntri))

    ! first compute the the smallest size of each triangle
    ! and store in diam_tri
    do count1=1,ntri
      Pa(:)=points(:,triang(1,count1))
      Pb(:)=points(:,triang(2,count1))
      Pc(:)=points(:,triang(3,count1))
      d_vect(1)=norm2(Pa-Pb)
      d_vect(2)=norm2(Pb-Pc)
      d_vect(3)=norm2(Pa-Pc)      
      diam_tri(count1)=minval(d_vect)**2
    enddo

    ! initialize flags that will identify the removed (redundant) points
    do count1=1,npoints
      flags(count1)=1
    enddo

    ! initialize tri_map, to locate all the triangles that contain certain point
    do count1=1,npoints
      do count2=1,nspace
        tri_map(count2,count1)=0
      enddo
    enddo
    ! fill tri_map
    do count1=1,ntri
      call add_tri(count1,triang(1,count1),tri_map,nspace,npoints)
      call add_tri(count1,triang(2,count1),tri_map,nspace,npoints)
      call add_tri(count1,triang(3,count1),tri_map,nspace,npoints)
    enddo

    ! this will be replace by search on nearest neigbours instead of n**2 algorithm
    do count2=1,npoints
      do count1=count2+1,npoints
        if (flags(count1).eq.1) then
          Pa=points(:,count1)
          Pb=points(:,count2)
          d=(Pa(1)-Pb(1))**2+(Pa(2)-Pb(2))**2+(Pa(3)-Pb(3))**2
          tol=eps*diam_tri(tri_map(1,count2))  !! aqui hay un error!!
          if (d.le.tol) then 
            flags(count1)=0
            call collapse_tri(count2,count1,tri_map,nspace,npoints)
          endif
        endif 
      enddo
    enddo
    npoints_out=sum(flags)

    ! now reconstruct the new triangles and points
    ! initialize triang_out
    do count1=1,ntri
      triang_out(1,count1)=0
      triang_out(2,count1)=0
      triang_out(3,count1)=0
    enddo

    icount=1
    do count1=1,npoints
      if (flags(count1).eq.1) then
        points_out(:,icount)=points(:,count1)
        icount=icount+1
        do count2=1,nspace
          if (tri_map(count2,count1).eq.0) then 
            exit
          else
            do count3=1,3
              if (triang_out(count3,tri_map(count2,count1)).eq.0) then
                triang_out(count3,tri_map(count2,count1))=icount-1
                exit
              endif
            enddo
          endif
        enddo
      endif
    enddo

    ! now we have to fix the orientation of the new triangles
    ! preserving the previous orientation
    do count1=1,ntri

      Pa1=points(:,triang(1,count1))
      Pb1=points(:,triang(2,count1))
      Pc1=points(:,triang(3,count1))
      call my_cross_v3(Pb1-Pa1, Pc1-Pa1, n1)
      Pa2=points_out(:,triang_out(1,count1))
      Pb2=points_out(:,triang_out(2,count1))
      Pc2=points_out(:,triang_out(3,count1))

      call my_cross_v3(Pb2-Pa2, Pc2-Pa2, n2)
!      write (*,*) 'dotprod: ',DOT_PRODUCT(n1,n2),n1,n2,Pb1-Pa1, Pc1-Pa1
!      read (*,*)

      if (DOT_PRODUCT(n1,n2)<0) then
        n_aux=triang_out(1,count1)
        triang_out(1,count1)=triang_out(2,count1)
        triang_out(2,count1)=n_aux
      endif
    enddo


return
end subroutine fix_tri_geom






subroutine my_cross_v3(a, b, c)
implicit none

    real ( kind = 8 ), intent(in) :: a(3),b(3)
    real ( kind = 8 ), intent(out) :: c(3)

        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)

end subroutine my_cross_v3






subroutine collapse_tri(point_num1,point_num2,tri_map,nspace,npoints)
implicit none

  !List of calling arguments
  integer, intent(in) :: point_num1,point_num2,nspace,npoints
  integer, intent(inout) :: tri_map(nspace,npoints)

  !List of local variables
  integer count1,count2

  do count1=1,nspace
    if (tri_map(count1,point_num2).eq.0) then
      exit
    else
      call add_tri(tri_map(count1,point_num2),point_num1,tri_map,nspace,npoints)
    endif
  enddo

return
end subroutine collapse_tri






subroutine add_tri(tri_num,point_num,tri_map,nspace,npoints)
implicit none

  !List of calling arguments
  integer, intent(in) :: tri_num,point_num,nspace,npoints
  integer, intent(inout) :: tri_map(nspace,npoints)

  !List of local variables
  integer count1

  do count1=1,nspace
    if (tri_map(count1,point_num).eq.0) then
      tri_map(count1,point_num)=tri_num
      exit
    elseif (count1.eq.nspace) then
      write (*,*) 'ERROR!! please, increase the variable nspace and try again..'
      stop
    endif
  enddo

return
end subroutine add_tri





!
!
!  This file contains the input routines for opening .a.tri files

subroutine open_tri_mem(filename,ntri,npoints)
implicit none

  !List of calling arguments
  character (len=*), intent(in) :: filename
  integer, intent(out) :: npoints,ntri
	
  !List of local variables
  integer umio,count1,count2,flag,aux
  integer :: ierror,norder

    open(UNIT=18, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ierror)

    read(18,*) npoints
    read(18,*) ntri
		
    close (18)

return
end subroutine open_tri_mem


subroutine tri2go3(triang,points,ntri,npoints,srcvals,srccoefs,norder,ns,wts)
implicit none

  !List of calling arguments
  integer, intent(in) :: ntri,npoints,ns,norder
  integer, intent(in) :: triang(3,ntri)
  real ( kind = 8 ), intent(in) :: points(3,npoints)
  real ( kind = 8 ), intent(out) :: srcvals(12,ns),wts(ns)
  real ( kind = 8 ), intent(out) :: srccoefs(9,ns)

  !List of local variables
  integer count1,count2,count3,count4,t1,t2,caux,icount
  integer npols,i1,i2,v1i,v2i,v3i
  real ( kind = 8 ) Pa(3),Pb(3),Pc(3),n_vect(3),aux_real,aux_vect(3)
  real ( kind = 8 ), allocatable :: rx(:),ry(:),rz(:)
  real ( kind = 8 ), allocatable :: uv(:,:),umatr(:,:),vmatr(:,:),w(:)
  real ( kind = 8 ), allocatable :: h_points(:),h_coefs(:)

  npols = (norder+1)*(norder+2)/2
  allocate(rx(npols),ry(npols),rz(npols))
  allocate(uv(2,npols),umatr(npols,npols),vmatr(npols,npols),w(npols))

  call vioreanu_simplex_quad(norder,npols,uv,umatr,vmatr,w)

  do count1=1,ntri
    Pa(:)=points(:,triang(1,count1))
    Pb(:)=points(:,triang(2,count1))
    Pc(:)=points(:,triang(3,count1))
    rx=(Pb(1)-Pa(1))*uv(1,:)+(Pc(1)-Pa(1))*uv(2,:)+Pa(1)
    ry=(Pb(2)-Pa(2))*uv(1,:)+(Pc(2)-Pa(2))*uv(2,:)+Pa(2)
    rz=(Pb(3)-Pa(3))*uv(1,:)+(Pc(3)-Pa(3))*uv(2,:)+Pa(3)
    i1=npols*(count1-1)+1
    i2=npols*count1
    srcvals(1,i1:i2)=rx
    srcvals(2,i1:i2)=ry
    srcvals(3,i1:i2)=rz
    do count2=1,npols
      srcvals(4,npols*(count1-1)+count2)=(Pb(1)-Pa(1))
      srcvals(5,npols*(count1-1)+count2)=(Pb(2)-Pa(2))
      srcvals(6,npols*(count1-1)+count2)=(Pb(3)-Pa(3))
      srcvals(7,npols*(count1-1)+count2)=(Pc(1)-Pa(1))
      srcvals(8,npols*(count1-1)+count2)=(Pc(2)-Pa(2))
      srcvals(9,npols*(count1-1)+count2)=(Pc(3)-Pa(3))

      call my_cross_v3(Pb-Pa, Pc-Pa, n_vect)
      n_vect=n_vect/norm2(n_vect)

      srcvals(10,npols*(count1-1)+count2)=n_vect(1)
      srcvals(11,npols*(count1-1)+count2)=n_vect(2)
      srcvals(12,npols*(count1-1)+count2)=n_vect(3)
    enddo
  enddo

	icount=1
	
  do count1=1,ntri
    do count2=1,npols
	  call my_cross_v3(srcvals(4:6,icount),srcvals(7:9,icount),aux_vect)
	    aux_real=sqrt(aux_vect(1)**2+aux_vect(2)**2+aux_vect(3)**2)
		wts(icount)=w(count2)*aux_real
		icount=icount+1
	enddo		
  enddo
  
  
  allocate(h_points(npols))
  allocate(h_coefs(npols))

  do count1=1,ntri
    h_points(:)=srcvals(1,(count1-1)*npols+1:count1*npols)
    call dmatvec(npols,npols,umatr,h_points,h_coefs)		
    srccoefs(1,(count1-1)*npols+1:count1*npols)=h_coefs(:)					
			
    h_points(:)=srcvals(2,(count1-1)*npols+1:count1*npols)
    call dmatvec(npols,npols,umatr,h_points,h_coefs)		
    srccoefs(2,(count1-1)*npols+1:count1*npols)=h_coefs(:)		

    h_points(:)=srcvals(3,(count1-1)*npols+1:count1*npols)
    call dmatvec(npols,npols,umatr,h_points,h_coefs)		
    srccoefs(3,(count1-1)*npols+1:count1*npols)=h_coefs(:)		

    h_points(:)=srcvals(4,(count1-1)*npols+1:count1*npols)
    call dmatvec(npols,npols,umatr,h_points,h_coefs)		
    srccoefs(4,(count1-1)*npols+1:count1*npols)=h_coefs(:)					
			
    h_points(:)=srcvals(5,(count1-1)*npols+1:count1*npols)
    call dmatvec(npols,npols,umatr,h_points,h_coefs)		
    srccoefs(5,(count1-1)*npols+1:count1*npols)=h_coefs(:)		

    h_points(:)=srcvals(6,(count1-1)*npols+1:count1*npols)
    call dmatvec(npols,npols,umatr,h_points,h_coefs)		
    srccoefs(6,(count1-1)*npols+1:count1*npols)=h_coefs(:)		

    h_points(:)=srcvals(7,(count1-1)*npols+1:count1*npols)
    call dmatvec(npols,npols,umatr,h_points,h_coefs)		
    srccoefs(7,(count1-1)*npols+1:count1*npols)=h_coefs(:)					
			
    h_points(:)=srcvals(8,(count1-1)*npols+1:count1*npols)
    call dmatvec(npols,npols,umatr,h_points,h_coefs)		
    srccoefs(8,(count1-1)*npols+1:count1*npols)=h_coefs(:)		

    h_points(:)=srcvals(9,(count1-1)*npols+1:count1*npols)
    call dmatvec(npols,npols,umatr,h_points,h_coefs)		
    srccoefs(9,(count1-1)*npols+1:count1*npols)=h_coefs(:)		
  enddo

return
end subroutine tri2go3



subroutine orthonormalize_all_v2(du,normal,ru,rv,ns)
implicit none

    !List of calling arguments
	integer, intent(in) :: ns
	real ( kind = 8 ), intent(in) :: du(3,ns), normal(3,ns)
	real ( kind = 8 ), intent(out) :: ru(3,ns), rv(3,ns)

	!List of local variables
	real ( kind = 8 ) aux
	integer count1

	do count1=1,ns
		call orthonormalize_v2(du(:,count1),normal(:,count1),ru(:,count1),rv(:,count1))
	enddo

return
end subroutine orthonormalize_all_v2

subroutine orthonormalize_v2(du,normal,ru,rv)
implicit none

    !List of calling arguments
	real ( kind = 8 ), intent(in) :: du(3), normal(3)
	real ( kind = 8 ), intent(out) :: ru(3), rv(3)

	!List of local variables
	real ( kind = 8 ) aux

	aux=sqrt(du(1)**2+du(2)**2+du(3)**2)
	
	ru(1)=du(1)/aux
	ru(2)=du(2)/aux
	ru(3)=du(3)/aux
		
	call my_cross_v3(normal, ru, rv)

return
end subroutine orthonormalize_v2








subroutine go32tri(srcvals,srccoefs,norder,ns,wts,N,npoints_est,triang,points,ntri,npoints)
implicit none

  !List of calling arguments
  real ( kind = 8 ), intent(in) :: srcvals(12,ns),wts(ns),srccoefs(9,ns)
  integer, intent(in) :: ns,norder,N,npoints_est
  integer, intent(out) :: npoints
  integer, intent(in) :: ntri
  integer, intent(out) :: triang(3,ntri)
  real ( kind = 8 ), intent(out) :: points(3,npoints_est)

  !List of local variables
  integer ntri_srcvals,npols
  integer n_vert_flat,n_tri_flat	!number of vertices and triangles of the resulting flat surface
  real ( kind = 8 ), allocatable :: vert_flat(:,:),u_array(:,:),v_array(:,:),srccoefs_aux(:)
  integer, allocatable :: tri_flat(:,:)
  real ( kind = 8 ), allocatable :: Px(:,:),Py(:,:),Pz(:,:)
  real ( kind = 8 ) uv(2),du,dv
  integer count0,count1,count2,icount1,icount2,icount3,n_aux

  npols = (norder+1)*(norder+2)/2

 ! (N input variable ) TOTAL NUMBER OF FLAT TRIANGLES IN THE PLOT, YOU MIGHT WANT TO DECREASE IT FOR VERY COMPLICATED GEOMETRIES
	
	write (*,*) 'total number of flat triangles per smooth triangle: ', N**2,'HOLAA'

!	n_aux=(N+1)*(N+2)/2
	allocate(u_array(N+1,N+1),v_array(N+1,N+1))
  write (*,*) 'medio'
	allocate(Px(N+1,N+1),Py(N+1,N+1),Pz(N+1,N+1))
write (*,*) 'ADIOSS'
	do count1=0,N
		do count2=0,N
			u_array(count1+1,count2+1)=real(count1,8)/real(N,8)
			v_array(count1+1,count2+1)=real(count2,8)/real(N,8)
!			write (*,*) u_array(count1+1,count2+1),v_array(count1+1,count2+1)
		enddo
	enddo
	du=1.0d0/(real(N,8))
	dv=1.0d0/(real(N,8))
	write (*,*) 'Dentro1',ns/npols
!	n_vert_flat=(N+1)**2*ntri
  ntri_srcvals=ns/npols
	n_tri_flat=N**2	*ntri_srcvals
	n_vert_flat=n_tri_flat*3
	allocate(vert_flat(3,n_vert_flat))
	allocate(tri_flat(3,n_tri_flat))
	allocate(srccoefs_aux(npols))
	write (*,*) 'Dentro2'

	icount1=1
	icount2=1
	icount3=1
	do count0=1,ntri_srcvals
	!here do the matvec to find the points
	  do count1=0,N
		do count2=0,N-count1
!			u_array=(/real(count1,8)/N,(real(count1,8)+1)/N,real(count1,8)/N,(real(count1,8)+1)/N /)
!			v_array=(/real(count1,8)/N,real(count1,8)/N,(real(count1,8)+1)/N,(real(count1,8)+1)/N /)
			uv(1)=u_array(count1+1,count2+1)
			uv(2)=v_array(count1+1,count2+1)
			srccoefs_aux(:)=srccoefs(1,(count0-1)*npols+1:count0*npols)
			call  koorn_evalexp(norder, npols, uv, srccoefs_aux, Px(count1+1,count2+1))
			srccoefs_aux(:)=srccoefs(2,(count0-1)*npols+1:count0*npols)
			call  koorn_evalexp(norder, npols, uv, srccoefs_aux, Py(count1+1,count2+1))
			srccoefs_aux(:)=srccoefs(3,(count0-1)*npols+1:count0*npols)
			call  koorn_evalexp(norder, npols, uv, srccoefs_aux, Pz(count1+1,count2+1))
		enddo
	  enddo

	  do count1=0,N-1
		do count2=0,N-1-count1
			vert_flat(1,icount1)=Px(count1+1,count2+1)
			vert_flat(2,icount1)=Py(count1+1,count2+1)
			vert_flat(3,icount1)=Pz(count1+1,count2+1)
			vert_flat(1,icount1+1)=Px(count1+1+1,count2+1)
			vert_flat(2,icount1+1)=Py(count1+1+1,count2+1)
			vert_flat(3,icount1+1)=Pz(count1+1+1,count2+1)		
			vert_flat(1,icount1+2)=Px(count1+1,count2+1+1)
			vert_flat(2,icount1+2)=Py(count1+1,count2+1+1)
			vert_flat(3,icount1+2)=Pz(count1+1,count2+1+1)
!			write (*,*) icount2
			tri_flat(1,icount2)=icount1
			tri_flat(2,icount2)=icount1+1
			tri_flat(3,icount2)=icount1+2			
			icount1=icount1+3
			icount2=icount2+1
		enddo
	  enddo
	  
	  do count1=1,N-1
		do count2=1,N-1-count1+1
			vert_flat(1,icount1)=Px(count1+1,count2+1)
			vert_flat(2,icount1)=Py(count1+1,count2+1)
			vert_flat(3,icount1)=Pz(count1+1,count2+1)
			vert_flat(1,icount1+1)=Px(count1+1-1,count2+1)
			vert_flat(2,icount1+1)=Py(count1+1-1,count2+1)
			vert_flat(3,icount1+1)=Pz(count1+1-1,count2+1)
			vert_flat(1,icount1+2)=Px(count1+1,count2+1-1)
			vert_flat(2,icount1+2)=Py(count1+1,count2+1-1)
			vert_flat(3,icount1+2)=Pz(count1+1,count2+1-1)		
!			write (*,*) icount2
			tri_flat(1,icount2)=icount1
			tri_flat(2,icount2)=icount1+1
			tri_flat(3,icount2)=icount1+2			
			icount1=icount1+3
			icount2=icount2+1
		enddo
	  enddo

	enddo

!do count1=1,n_vert_flat
!  write (*,*) vert_flat(:,count1)
!enddo
!do count1=1,n_tri_flat
!  write (*,*) tri_flat(:,count1)
!enddo
!write (*,*) 'no limpio'
!stop
write (*,*) 'almost done'
  call fix_tri_geom(n_tri_flat,n_vert_flat,vert_flat,tri_flat,npoints,points,triang)
write (*,*) 'done'
return
end subroutine go32tri







subroutine open_go3_4_rwg(fname,npts_aux,npatches_aux,norder,npatches, &
    norders,ixyzs,&
    iptype,npts,srcvals,srccoefs,wts,N,npoints_est,nrwg,rwg,vert_flat,&
    tri_flat,n_vert_flat,n_tri_flat)
  implicit none
  !
  ! This routine ...
  !

  ! List of calling arguments
    character (len=*), intent(in) :: fname
    integer, intent(in) :: npts_aux,npatches_aux,norder,npatches
    integer, intent(out) :: norders(npatches)
    integer, intent(out) :: ixyzs(npatches+1)
    integer, intent(out) :: iptype(npatches)
    integer, intent(in) :: npts
    real ( kind = 8 ), intent(out) :: srcvals(12,npts)
    real ( kind = 8 ), intent(out) :: srccoefs(9,npts)
    real ( kind = 8 ), intent(out) :: wts(npts)
    integer, intent(in) :: N,npoints_est
    integer, intent(out) :: nrwg
    integer, intent(out) :: rwg(4,3*npatches)
    real ( kind = 8 ), intent(out) :: vert_flat(3,npoints_est)
    integer, intent(out) :: tri_flat(3,npatches)
    integer, intent(out) :: n_vert_flat,n_tri_flat

    !List of local variables
    integer, allocatable :: norders_aux(:),ixyzs_aux(:),iptype_aux(:)
    real ( kind = 8 ), allocatable :: srcvals_aux(:,:),srccoefs_aux(:,:),wts_aux(:)
    integer npols,norder_aux,i

  allocate(srcvals_aux(12,npts_aux),srccoefs_aux(9,npts_aux))
  allocate(ixyzs_aux(npatches_aux+1),iptype_aux(npatches_aux),norders_aux(npatches_aux))
  allocate(wts_aux(npts_aux))

      call open_gov3_geometry(fname,npatches_aux,norders_aux,ixyzs_aux,&
       &iptype_aux,npts_aux,srcvals_aux,srccoefs_aux,wts_aux)
     
      norder_aux=norders_aux(1)

      n_tri_flat=npatches
	  
      call go32tri(srcvals_aux,srccoefs_aux,norder_aux,npts_aux,wts_aux,N,&
       &npoints_est,tri_flat,vert_flat,n_tri_flat,n_vert_flat)
      
       call rwg_basis(n_tri_flat,n_vert_flat,vert_flat,tri_flat,rwg,nrwg)

      call tri2go3(tri_flat,vert_flat,n_tri_flat,n_vert_flat,srcvals,srccoefs,norder,npts,wts)

      npols = (norder+1)*(norder+2)/2
      do i = 1,npatches
        norders(i) = norder
        iptype(i) = 1
        ixyzs(i) = (i-1)*npols + 1
      enddo
      ixyzs(npatches+1) = npts+1

return
end subroutine open_go3_4_rwg




subroutine open_go3_4_rwg_v2(srcvals_aux,srccoefs_aux,wts_aux,ixyzs_aux,iptype_aux,norders_aux,&
     & npts_aux,npatches_aux,norder,npatches,norders,ixyzs,&
     &iptype,npts,srcvals,srccoefs,wts,N,npoints_est,nrwg,rwg,vert_flat,&
     &tri_flat,n_vert_flat,n_tri_flat)
implicit none

    !List of calling arguments
    integer, intent(in) :: npts_aux,npatches_aux,norder,npatches
    integer, intent(out) :: norders(npatches)
    integer, intent(out) :: ixyzs(npatches+1)
    integer, intent(out) :: iptype(npatches)

    integer, intent(in) :: norders_aux(npatches_aux)
    integer, intent(in) :: iptype_aux(npatches_aux)
    integer, intent(in) :: ixyzs_aux(npatches_aux+1)
    
    integer, intent(in) :: npts

    real ( kind = 8 ), intent(out) :: srcvals(12,npts)
    real ( kind = 8 ), intent(out) :: srccoefs(9,npts)
    real ( kind = 8 ), intent(out) :: wts(npts)

    real ( kind = 8 ), intent(in) :: srcvals_aux(12,npts_aux)
    real ( kind = 8 ), intent(in) :: srccoefs_aux(9,npts_aux)
    real ( kind = 8 ), intent(in) :: wts_aux(npts_aux)


    integer, intent(in) :: N,npoints_est
    integer, intent(out) :: nrwg
    integer, intent(out) :: rwg(4,3*npatches)
    real ( kind = 8 ), intent(out) :: vert_flat(3,npoints_est)
    integer, intent(out) :: tri_flat(3,npatches)
    integer, intent(out) :: n_vert_flat,n_tri_flat

    !List of local variables
    integer npols,norder_aux,i

!      allocate(srcvals_aux(12,npts_aux),srccoefs_aux(9,npts_aux))
!      allocate(ixyzs_aux(npatches_aux+1),iptype_aux(npatches_aux),norders_aux(npatches_aux))
!      allocate(wts_aux(npts_aux))

!      call open_gov3_geometry(fname,npatches_aux,norders_aux,ixyzs_aux,&
!       &iptype_aux,npts_aux,srcvals_aux,srccoefs_aux,wts_aux)
     
      norder_aux=norders_aux(1)

      n_tri_flat=npatches
	  
      call go32tri(srcvals_aux,srccoefs_aux,norder_aux,npts_aux,wts_aux,N,&
       &npoints_est,tri_flat,vert_flat,n_tri_flat,n_vert_flat)
      call rwg_basis(n_tri_flat,n_vert_flat,vert_flat,tri_flat,rwg,nrwg)

      call tri2go3(tri_flat,vert_flat,n_tri_flat,n_vert_flat,srcvals,srccoefs,norder,npts,wts)

      npols = (norder+1)*(norder+2)/2
      do i = 1,npatches
        norders(i) = norder
        iptype(i) = 1
        ixyzs(i) = (i-1)*npols + 1
      enddo
      ixyzs(npatches+1) = npts+1

return
end subroutine open_go3_4_rwg_v2
