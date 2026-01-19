

subroutine cisurf_quad2flat(Geometry1)
  use ModType_Smooth_Surface
  implicit double precision (a-h,o-z)

  type (Geometry) :: Geometry1
  
  integer *8 :: js(10)
  double precision :: xyz1(10), xyz2(10), xyz3(10)
  
  !
  ! this subroutine flattens the quadratic triangles that were loaded
  ! into 1st order, flat triangles
  !

  ntri = Geometry1%ntri_sk

  do i = 1,ntri

    do j = 1,6
      js(j) = Geometry1%Tri(j,i)
    end do

    do k = 1,3
      xyz1(k) = Geometry1%Points(k,js(1))
      xyz2(k) = Geometry1%Points(k,js(2))
      xyz3(k) = Geometry1%Points(k,js(3))
    end do

    do k = 1,3
      Geometry1%Points(k,js(4)) = (xyz1(k) + xyz2(k))/2
      Geometry1%Points(k,js(5)) = (xyz2(k) + xyz3(k))/2
      Geometry1%Points(k,js(6)) = (xyz3(k) + xyz1(k))/2
    end do
    
  end do

  Geometry1%ifflat = 1

  return
end subroutine cisurf_quad2flat





subroutine funcion_skeleton(Geometry1)
  use ModType_Smooth_Surface
  implicit none

  !List of calling arguments
  integer *8 :: norder_skel
  type (Geometry), intent(inout) :: Geometry1
  integer *8 :: nsk


  !List of local variables
  double precision, allocatable :: UV(:,:),w(:),U(:),V(:)
  double precision :: P1(3),P2(3),P3(3),P4(3),P5(3),P6(3)
  double precision, allocatable :: F_x(:), F_y(:), F_z(:), dS(:)
  double precision, allocatable :: nP_x(:), nP_y(:), nP_z(:)
  integer *8 count, itype, npols

  double precision, allocatable :: umatr(:,:), vmatr(:,:)

  norder_skel = Geometry1%norder_skel  
  nsk = (norder_skel+1)*(norder_skel+2)/2

  allocate(umatr(nsk,nsk),vmatr(nsk,nsk))
  allocate( UV(2,nsk), w(nsk) )
  allocate( F_x(nsk), F_y(nsk), F_z(nsk), dS(nsk) )
  allocate( nP_x(nsk), nP_y(nsk), nP_z(nsk) )
  
  ! get the specified quadrature nodes and weights on the simplex of
  ! order norder_skel
  
  npols = nsk 
  call vioreanu_simplex_quad(norder_skel,npols, UV, umatr,vmatr,w)
  U = UV(1,:)
  V = UV(2,:)
   

  
  if (allocated(Geometry1%skeleton_Points)) then
    deallocate(Geometry1%skeleton_Points)
  endif

  if (allocated(Geometry1%skeleton_w)) then
    deallocate(Geometry1%skeleton_w)
  endif

  if (allocated(Geometry1%skeleton_N)) then
    deallocate(Geometry1%skeleton_N)
  endif

  
  allocate(Geometry1%skeleton_Points(3,Geometry1%n_Sk_points))
  allocate(Geometry1%skeleton_w(Geometry1%n_Sk_points))
  allocate(Geometry1%skeleton_N(3,Geometry1%n_Sk_points))

  do count=1,Geometry1%ntri_sk

    P1=Geometry1%Points(:,Geometry1%Tri(1,count))
    P2=Geometry1%Points(:,Geometry1%Tri(2,count))
    P3=Geometry1%Points(:,Geometry1%Tri(3,count))
    P4=Geometry1%Points(:,Geometry1%Tri(4,count))
    P5=Geometry1%Points(:,Geometry1%Tri(5,count))
    P6=Geometry1%Points(:,Geometry1%Tri(6,count))

    call eval_quadratic_patch(P1,P2,P3,P4,P5,P6,U,V,F_x,F_y,F_z,nP_x,nP_y,nP_z,dS,nsk)

    Geometry1%skeleton_Points(1,(count-1)*nsk+1:(count)*nsk)=F_x
    Geometry1%skeleton_Points(2,(count-1)*nsk+1:(count)*nsk)=F_y
    Geometry1%skeleton_Points(3,(count-1)*nsk+1:(count)*nsk)=F_z
    Geometry1%skeleton_w((count-1)*nsk+1:(count)*nsk)=w*dS

    Geometry1%skeleton_N(1,(count-1)*nsk+1:(count)*nsk)=nP_x
    Geometry1%skeleton_N(2,(count-1)*nsk+1:(count)*nsk)=nP_y
    Geometry1%skeleton_N(3,(count-1)*nsk+1:(count)*nsk)=nP_z
  enddo

  return
end subroutine funcion_skeleton






subroutine eval_quadratic_patch(P1,P2,P3,P4,P5,P6,U,V, &
    F_x, F_y, F_z, nP_x, nP_y, nP_z, dS, n_order)
  implicit none

  !
  ! this routine evaluates 
  !

  !List of calling arguments
  integer *8, intent(in) :: n_order
  double precision, intent(in) :: P1(3),P2(3),P3(3),P4(3),P5(3),P6(3)
  double precision, intent(in) :: U(n_order),V(n_order)
  double precision, intent(out) :: F_x(n_order),F_y(n_order),F_z(n_order),dS(n_order)
  double precision, intent(out) :: nP_x(n_order),nP_y(n_order),nP_z(n_order)

  !List of local variables
  double precision coef_x(6),coef_y(6),coef_z(6),U_x,U_y,U_z,V_x,V_y,V_z
  integer *8 count

  coef_x(1)=P1(1)
  coef_x(2)=-3*P1(1)-P2(1)+4*P4(1)
  coef_x(3)=-3*P1(1)-P3(1)+4*P6(1)
  coef_x(4)=2*P1(1)+2*P2(1)-4*P4(1)
  coef_x(5)=2*P1(1)+2*P3(1)-4*P6(1)
  coef_x(6)=4*P1(1)-4*P4(1)+4*P5(1)-4*P6(1)

  coef_y(1)=P1(2)
  coef_y(2)=-3*P1(2)-P2(2)+4*P4(2)
  coef_y(3)=-3*P1(2)-P3(2)+4*P6(2)
  coef_y(4)=2*P1(2)+2*P2(2)-4*P4(2)
  coef_y(5)=2*P1(2)+2*P3(2)-4*P6(2)
  coef_y(6)=4*P1(2)-4*P4(2)+4*P5(2)-4*P6(2)

  coef_z(1)=P1(3)
  coef_z(2)=-3*P1(3)-P2(3)+4*P4(3)
  coef_z(3)=-3*P1(3)-P3(3)+4*P6(3)
  coef_z(4)=2*P1(3)+2*P2(3)-4*P4(3)
  coef_z(5)=2*P1(3)+2*P3(3)-4*P6(3)
  coef_z(6)=4*P1(3)-4*P4(3)+4*P5(3)-4*P6(3)

  do count=1,n_order
    F_x(count)=coef_x(1)+coef_x(2)*U(count)+coef_x(3)*V(count)+coef_x(4)&
        &*U(count)**2+coef_x(5)*V(count)**2+coef_x(6)*U(count)*V(count)
    F_y(count)=coef_y(1)+coef_y(2)*U(count)+coef_y(3)*V(count)+coef_y(4)&
        &*U(count)**2+coef_y(5)*V(count)**2+coef_y(6)*U(count)*V(count)
    F_z(count)=coef_z(1)+coef_z(2)*U(count)+coef_z(3)*V(count)+coef_z(4)&
        &*U(count)**2+coef_z(5)*V(count)**2+coef_z(6)*U(count)*V(count)
    U_x=coef_x(2)+2*coef_x(4)*U(count)+coef_x(6)*V(count)
    U_y=coef_y(2)+2*coef_y(4)*U(count)+coef_y(6)*V(count)
    U_z=coef_z(2)+2*coef_z(4)*U(count)+coef_z(6)*V(count)
    V_x=coef_x(3)+2*coef_x(5)*V(count)+coef_x(6)*U(count)
    V_y=coef_y(3)+2*coef_y(5)*V(count)+coef_y(6)*U(count)
    V_z=coef_z(3)+2*coef_z(5)*V(count)+coef_z(6)*U(count)
    nP_x(count)=U_y*V_z-U_z*V_y;
    nP_y(count)=U_z*V_x-U_x*V_z;
    nP_z(count)=U_x*V_y-U_y*V_x;
    dS(count)=sqrt(nP_x(count)**2+nP_y(count)**2+nP_z(count)**2);
    nP_x(count)=nP_x(count)/dS(count)
    nP_y(count)=nP_y(count)/dS(count)
    nP_z(count)=nP_z(count)/dS(count)
  enddo
  return
end subroutine eval_quadratic_patch






subroutine funcion_normal_vert(Geometry1)
  use ModType_Smooth_Surface
  implicit none

  !List of calling arguments
  type (Geometry), intent(inout) :: Geometry1

  double precision find_angle

  
  !List of local variables
  type (My_cell) My_cell1
  integer *8 count, n_order
  double precision P1(3),P2(3),P3(3),P4(3),P5(3),P6(3)
  double precision U(3),V(3)
  double precision F_x(3),F_y(3),F_z(3),dS(3)
  double precision U_x(3),U_y(3),U_z(3),V_x(3),V_y(3),V_z(3)
  double precision angle_vertex(3)
  double precision nP_x(3),nP_y(3),nP_z(3)
  double precision current_Normal_vert(3)

  U= (/0.0d0, 1.0d0, 0.0d0/)
  V= (/0.0d0, 0.0d0, 1.0d0/)
  n_order=3
  My_cell1%n_Cell=Geometry1%npoints
  allocate(My_cell1%Var_Mat(Geometry1%npoints))
  do count=1,Geometry1%npoints
    My_cell1%Var_Mat(count)%n_Mat=0
    My_cell1%Var_Mat(count)%current_n_Mat=1
  enddo
  do count=1,Geometry1%ntri
    My_cell1%Var_Mat(Geometry1%Tri(1,count))%n_Mat=My_cell1%Var_Mat(Geometry1%Tri(1,count))%n_Mat+1
    My_cell1%Var_Mat(Geometry1%Tri(2,count))%n_Mat=My_cell1%Var_Mat(Geometry1%Tri(2,count))%n_Mat+1
    My_cell1%Var_Mat(Geometry1%Tri(3,count))%n_Mat=My_cell1%Var_Mat(Geometry1%Tri(3,count))%n_Mat+1
  enddo
  do count=1,Geometry1%npoints
    if (My_cell1%Var_Mat(count)%n_Mat>0) then
      allocate(My_cell1%Var_Mat(count)%Mat(4,My_cell1%Var_Mat(count)%n_Mat))
    endif
  enddo
  do count=1,Geometry1%ntri
    P1=Geometry1%Points(:,Geometry1%Tri(1,count))
    P2=Geometry1%Points(:,Geometry1%Tri(2,count))
    P3=Geometry1%Points(:,Geometry1%Tri(3,count))
    P4=Geometry1%Points(:,Geometry1%Tri(4,count))
    P5=Geometry1%Points(:,Geometry1%Tri(5,count))
    P6=Geometry1%Points(:,Geometry1%Tri(6,count))
    !            call eval_quadratic_patch(P1,P2,P3,P4,P5,P6,U,V,F_x,F_y,F_z,nP_x,nP_y,nP_z,dS,n_order)

    call eval_quadratic_patch_UV(P1,P2,P3,P4,P5,P6,U,V,F_x,F_y,F_z,U_x,U_y,U_z,V_x,V_y,V_z,nP_x,nP_y,nP_z,dS,n_order)


    My_cell1%Var_Mat(Geometry1%Tri(1,count))%Mat(1:3,My_cell1%Var_Mat(Geometry1%Tri(1,count))%current_n_Mat)&
        &=(/nP_x(1), nP_y(1), nP_z(1)/)
    angle_vertex(1)=find_angle(U_x(1),U_y(1),U_z(1),V_x(1),V_y(1),V_z(1))
    My_cell1%Var_Mat(Geometry1%Tri(1,count))%Mat(4,My_cell1%Var_Mat&
        &(Geometry1%Tri(1,count))%current_n_Mat)=angle_vertex(1)
    My_cell1%Var_Mat(Geometry1%Tri(1,count))%current_n_Mat=My_cell1%Var_Mat&
        &(Geometry1%Tri(1,count))%current_n_Mat+1

    !            write (*,*) Geometry1%Points(:,Geometry1%Tri(1,count))
    !            write (*,*) U_x(1),U_y(1),U_z(1)
    !            write (*,*) V_x(1),V_y(1),V_z(1)
    !            write (*,*) angle_vertex(1)/3.14159268d0

    My_cell1%Var_Mat(Geometry1%Tri(2,count))%Mat(1:3,My_cell1%Var_Mat&
        &(Geometry1%Tri(2,count))%current_n_Mat)=(/nP_x(2), nP_y(2), nP_z(2)/)
    angle_vertex(2)=find_angle(-U_x(2),-U_y(2),-U_z(2),V_x(2)-U_x(2),V_y(2)-U_y(2),V_z(2)-U_z(2))
    My_cell1%Var_Mat(Geometry1%Tri(2,count))%Mat(4,My_cell1%Var_Mat&
        &(Geometry1%Tri(2,count))%current_n_Mat)=angle_vertex(2)
    My_cell1%Var_Mat(Geometry1%Tri(2,count))%current_n_Mat=My_cell1%Var_Mat&
        &(Geometry1%Tri(2,count))%current_n_Mat+1

    !            write (*,*) Geometry1%Points(:,Geometry1%Tri(2,count))
    !            write (*,*) -U_x(2),-U_y(2),-U_z(2)
    !            write (*,*) V_x(2)-U_x(2),V_y(2)-U_y(2),V_z(2)-U_z(2)
    !            write (*,*) angle_vertex(2)/3.14159268d0


    My_cell1%Var_Mat(Geometry1%Tri(3,count))%Mat(1:3,My_cell1%Var_Mat&
        &(Geometry1%Tri(3,count))%current_n_Mat)=(/nP_x(3), nP_y(3), nP_z(3)/)
    angle_vertex(3)=find_angle(-V_x(3)+U_x(3),-V_y(3)+U_y(3),-V_z(3)+U_z(3),-V_x(3),-V_y(3),-V_z(3))
    My_cell1%Var_Mat(Geometry1%Tri(3,count))%Mat(4,My_cell1%Var_Mat&
        &(Geometry1%Tri(3,count))%current_n_Mat)=angle_vertex(3)
    My_cell1%Var_Mat(Geometry1%Tri(3,count))%current_n_Mat=My_cell1%Var_Mat&
        &(Geometry1%Tri(3,count))%current_n_Mat+1

    !            write (*,*) Geometry1%Points(:,Geometry1%Tri(3,count))
    !            write (*,*) -V_x(3)+U_x(3),-V_y(3)+U_y(3),-V_z(3)+U_z(3)
    !            write (*,*) -V_x(3),-V_y(3),-V_z(3)
    !            write (*,*) angle_vertex(3)/3.14159268d0

    !            read (*,*)

  enddo
  if (allocated(Geometry1%Normal_Vert)) then
    deallocate(Geometry1%Normal_Vert)
  endif
  allocate(Geometry1%Normal_Vert(3,Geometry1%npoints))
  do count=1,Geometry1%npoints
    if (My_cell1%Var_Mat(count)%n_Mat>0) then
      !                write (*,*) count
      !                write (*,*) My_cell1%Var_Mat(count)%Mat(1:3,:)
      !                write (*,*) count
      !                write (*,*) My_cell1%Var_Mat(count)%Mat(4,:)
      !                write (*,*) 'n point: ',count
      !                write (*,*) sum(My_cell1%Var_Mat(count)%Mat(4,:))/3.14159268d0
      !                write (*,*)
      !                write (*,*) Geometry1%Points(:,count)
      !                read (*,*)
      !                write (*,*) 'n_mat: ',My_cell1%Var_Mat(count)%n_Mat
      !                write (*,*) My_cell1%Var_Mat(count)%Mat(4,:)
      !                write (*,*) My_cell1%Var_Mat(count)%Mat(1:3,:)

      call mimean(My_cell1%Var_Mat(count)%Mat,My_cell1%Var_Mat(count)%n_Mat,current_Normal_vert)
      Geometry1%Normal_Vert(:,count)=current_Normal_vert
      !                write (*,*) 'aqui: ', current_Normal_vert
      !                read (*,*)
    else
      Geometry1%Normal_Vert(:,count)=(/0.d0, 0.d0, 0.d0/)
    endif
  enddo
  !        write (*,*) Geometry1%Normal_Vert

  return
end subroutine funcion_normal_vert





subroutine mimean(All_Normals,n_Normals,Current_Normal)
  implicit none

  !List of calling arguments
  integer *8, intent(in) :: n_Normals
  double precision, intent(in) :: All_Normals(4,n_Normals)
  double precision, intent(out) :: Current_Normal(3)

  !List of local variables
  integer *8 count
  double precision coef(n_Normals),sum_coef

  current_Normal=(/ 0.0d0, 0.0d0, 0.0d0 /)
  sum_coef=0
  do count=1,n_Normals
    sum_coef=sum_coef+All_Normals(4,count)
  enddo
  do count=1,n_Normals
    Current_Normal=Current_Normal+All_Normals(1:3,count)*All_Normals(4,count)
  enddo
  Current_Normal=Current_Normal/sum_coef
  !        write (*,*) 'all normals: ', All_Normals(1:3,:)
  !        write (*,*) 'all coefs: ', All_Normals(4,:)
  !        write (*,*) Current_Normal,sum_coef
  !        read (*,*)
  return
end subroutine mimean





  
subroutine eval_quadratic_patch_UV(P1,P2,P3,P4,P5,P6,U,V, &
    F_x,F_y,F_z,U_x,U_y,U_z,V_x,V_y,V_z,nP_x,nP_y,nP_z,dS,n_order)
    implicit none

    !
    ! this routine ...
    !
    
    !List of calling arguments
    integer *8, intent(in) :: n_order
    double precision, intent(in) :: P1(3),P2(3),P3(3),P4(3),P5(3),P6(3)
    double precision, intent(in) :: U(n_order),V(n_order)
    double precision, intent(out) :: F_x(n_order),F_y(n_order),F_z(n_order),dS(n_order)
    double precision, intent(out) :: nP_x(n_order),nP_y(n_order),nP_z(n_order)
    double precision, intent(out) :: U_x(n_order),U_y(n_order),U_z(n_order),V_x(n_order),V_y(n_order),V_z(n_order)

    !List of local variables
    double precision coef_x(6),coef_y(6),coef_z(6)
    integer *8 count

    coef_x(1)=P1(1)
    coef_x(2)=-3*P1(1)-P2(1)+4*P4(1)
    coef_x(3)=-3*P1(1)-P3(1)+4*P6(1)
    coef_x(4)=2*P1(1)+2*P2(1)-4*P4(1)
    coef_x(5)=2*P1(1)+2*P3(1)-4*P6(1)
    coef_x(6)=4*P1(1)-4*P4(1)+4*P5(1)-4*P6(1)

    coef_y(1)=P1(2)
    coef_y(2)=-3*P1(2)-P2(2)+4*P4(2)
    coef_y(3)=-3*P1(2)-P3(2)+4*P6(2)
    coef_y(4)=2*P1(2)+2*P2(2)-4*P4(2)
    coef_y(5)=2*P1(2)+2*P3(2)-4*P6(2)
    coef_y(6)=4*P1(2)-4*P4(2)+4*P5(2)-4*P6(2)

    coef_z(1)=P1(3)
    coef_z(2)=-3*P1(3)-P2(3)+4*P4(3)
    coef_z(3)=-3*P1(3)-P3(3)+4*P6(3)
    coef_z(4)=2*P1(3)+2*P2(3)-4*P4(3)
    coef_z(5)=2*P1(3)+2*P3(3)-4*P6(3)
    coef_z(6)=4*P1(3)-4*P4(3)+4*P5(3)-4*P6(3)

    do count=1,n_order
      F_x(count)=coef_x(1)+coef_x(2)*U(count)+coef_x(3)*V(count)+coef_x(4)&
          &*U(count)**2+coef_x(5)*V(count)**2+coef_x(6)*U(count)*V(count)
      F_y(count)=coef_y(1)+coef_y(2)*U(count)+coef_y(3)*V(count)+coef_y(4)&
          &*U(count)**2+coef_y(5)*V(count)**2+coef_y(6)*U(count)*V(count)
      F_z(count)=coef_z(1)+coef_z(2)*U(count)+coef_z(3)*V(count)+coef_z(4)&
          &*U(count)**2+coef_z(5)*V(count)**2+coef_z(6)*U(count)*V(count)
      U_x(count)=coef_x(2)+2*coef_x(4)*U(count)+coef_x(6)*V(count)
      U_y(count)=coef_y(2)+2*coef_y(4)*U(count)+coef_y(6)*V(count)
      U_z(count)=coef_z(2)+2*coef_z(4)*U(count)+coef_z(6)*V(count)
      V_x(count)=coef_x(3)+2*coef_x(5)*V(count)+coef_x(6)*U(count)
      V_y(count)=coef_y(3)+2*coef_y(5)*V(count)+coef_y(6)*U(count)
      V_z(count)=coef_z(3)+2*coef_z(5)*V(count)+coef_z(6)*U(count)
      nP_x(count)=U_y(count)*V_z(count)-U_z(count)*V_y(count);
      nP_y(count)=U_z(count)*V_x(count)-U_x(count)*V_z(count);
      nP_z(count)=U_x(count)*V_y(count)-U_y(count)*V_x(count);
      dS(count)=sqrt(nP_x(count)**2+nP_y(count)**2+nP_z(count)**2);
      nP_x(count)=nP_x(count)/dS(count)
      nP_y(count)=nP_y(count)/dS(count)
      nP_z(count)=nP_z(count)/dS(count)
    enddo
    return
  end subroutine eval_quadratic_patch_UV





  subroutine crossproduct(a, b, c)
    double precision :: a(3), b(3), c(3)
    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
    return
  end subroutine crossproduct
  




  function find_angle(u_x, u_y, u_z, v_x, v_y, v_z)
    double precision, intent(in) :: u_x, u_y, u_z, v_x, v_y, v_z
    double precision find_angle
    double precision num,denom
    num = u_x*v_x+u_y*v_y+u_z*v_z
    denom = sqrt(u_x**2+u_y**2+u_z**2)*sqrt(v_x**2+v_y**2+v_z**2)
    find_angle = acos(num/denom)
  end function find_angle
  
