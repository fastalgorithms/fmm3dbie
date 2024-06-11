Module Mod_Smooth_Surface
  
  use Mod_Feval
  use ModType_Smooth_Surface
  
  implicit none
  
contains
  
  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_cos(v,All_Normals,coef,n_Normals)
    implicit none

    !List of calling arguments
    integer, intent(in) :: n_Normals
    double precision, intent(out) :: coef
    double precision, intent(in) :: All_Normals(3,n_Normals), v(3)

    !List of local variables
    integer count
    double precision n_times, tol, my_dot, norm1, norm2
    n_times=0
    tol=1.0d-10
    norm1=sqrt(v(1)**2+v(2)**2+v(3)**2)
    do count=1,n_Normals
      norm2=sqrt(All_Normals(1,count)**2+All_Normals(2,count)**2+All_Normals(3,count)**2)
      my_dot=dot_product(v,All_Normals(:,count))/(norm1*norm2)
      if (my_dot<0) then
        my_dot=0
      endif
      n_times=n_times+my_dot**2
    enddo
    coef=1/n_times
    return
  end subroutine find_cos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine refine_geometry(Geometry1,n_order_sf)
    implicit none

    !List of calling arguments
    type (Geometry), intent(inout) :: Geometry1
    integer , intent(in) :: n_order_sf

    !List of local variables
    integer count,contador_indices
    double precision, allocatable :: Points(:,:), Normal_Vert(:,:)
    integer, allocatable :: Tri(:,:)
    double precision P1(3),P2(3),P3(3),P4(3),P5(3),P6(3)
    double precision Pa(3),Pb(3),Pc(3),Pd(3),Pe(3),Pf(3),Pg(3),Ph(3),Pi(3)
    double precision Nor1(3),Nor2(3),Nor3(3),Nor4(3),Nor5(3),Nor6(3)
    double precision Nor_a(3),Nor_b(3),Nor_c(3),Nor_d(3),Nor_e(3),Nor_f(3),Nor_g(3),Nor_h(3),Nor_i(3)
    double precision U(9),V(9)
    double precision F_x(9),F_y(9),F_z(9),dS(9)
    double precision nP_x(9),nP_y(9),nP_z(9)
    double precision U_x(9),U_y(9),U_z(9),V_x(9),V_y(9),V_z(9)
    integer m,N,n_order_aux

    allocate(Points(3,Geometry1%ntri*15))
    allocate(Normal_Vert(3,Geometry1%ntri*15))
    allocate(Tri(6,Geometry1%ntri*4))
    contador_indices=1
    n_order_aux=9
    U=[0.2500d0,0.7500d0,0d0,0.2500d0,0.5000d0,0.7500d0,0.2500d0,0d0,0.2500d0]
    V=[0d0,0d0,0.2500d0,0.2500d0,0.2500d0,0.2500d0,0.5000d0,0.7500d0,0.7500d0]
    do count=1,Geometry1%ntri
      P1=Geometry1%Points(:,Geometry1%Tri(1,count))
      P2=Geometry1%Points(:,Geometry1%Tri(2,count))
      P3=Geometry1%Points(:,Geometry1%Tri(3,count))
      P4=Geometry1%Points(:,Geometry1%Tri(4,count))
      P5=Geometry1%Points(:,Geometry1%Tri(5,count))
      P6=Geometry1%Points(:,Geometry1%Tri(6,count))
      call eval_quadratic_patch_UV(P1,P2,P3,P4,P5,P6,U,V,F_x,F_y,F_z,U_x,U_y,U_z,V_x,V_y,V_z,nP_x,nP_y,nP_z,dS,n_order_aux)
      Pa=[F_x(1),F_y(1),F_z(1)]
      Pb=[F_x(2),F_y(2),F_z(2)]
      Pc=[F_x(3),F_y(3),F_z(3)]
      Pd=[F_x(4),F_y(4),F_z(4)]
      Pe=[F_x(5),F_y(5),F_z(5)]
      Pf=[F_x(6),F_y(6),F_z(6)]
      Pg=[F_x(7),F_y(7),F_z(7)]
      Ph=[F_x(8),F_y(8),F_z(8)]
      Pi=[F_x(9),F_y(9),F_z(9)]
      Nor1=Geometry1%Normal_Vert(:,Geometry1%Tri(1,count))
      Nor2=Geometry1%Normal_Vert(:,Geometry1%Tri(2,count))
      Nor3=Geometry1%Normal_Vert(:,Geometry1%Tri(3,count))
      Nor4=(Nor1+Nor2)/2.0d0
      Nor5=(Nor2+Nor3)/2.0d0
      Nor6=(Nor3+Nor1)/2.0d0
      Nor_a=(Nor1+Nor4)/2.0d0
      Nor_b=(Nor4+Nor2)/2.0d0
      Nor_c=(Nor1+Nor6)/2.0d0
      Nor_d=(Nor4+Nor6)/2.0d0
      Nor_e=(Nor4+Nor5)/2.0d0
      Nor_f=(Nor5+Nor2)/2.0d0
      Nor_g=(Nor5+Nor6)/2.0d0
      Nor_h=(Nor3+Nor6)/2.0d0
      Nor_i=(Nor3+Nor5)/2.0d0
      Points(:,contador_indices)=P1
      Points(:,contador_indices+1)=Pa
      Points(:,contador_indices+2)=P4
      Points(:,contador_indices+3)=Pb
      Points(:,contador_indices+4)=P2
      Points(:,contador_indices+5)=Pc
      Points(:,contador_indices+6)=Pd
      Points(:,contador_indices+7)=Pe
      Points(:,contador_indices+8)=Pf
      Points(:,contador_indices+9)=P6
      Points(:,contador_indices+10)=Pg
      Points(:,contador_indices+11)=P5
      Points(:,contador_indices+12)=Ph
      Points(:,contador_indices+13)=Pi
      Points(:,contador_indices+14)=P3

      Normal_Vert(:,contador_indices)=Nor1
      Normal_Vert(:,contador_indices+1)=Nor_a
      Normal_Vert(:,contador_indices+2)=Nor4
      Normal_Vert(:,contador_indices+3)=Nor_b
      Normal_Vert(:,contador_indices+4)=Nor2
      Normal_Vert(:,contador_indices+5)=Nor_c
      Normal_Vert(:,contador_indices+6)=Nor_d
      Normal_Vert(:,contador_indices+7)=Nor_e
      Normal_Vert(:,contador_indices+8)=Nor_f
      Normal_Vert(:,contador_indices+9)=Nor6
      Normal_Vert(:,contador_indices+10)=Nor_g
      Normal_Vert(:,contador_indices+11)=Nor5
      Normal_Vert(:,contador_indices+12)=Nor_h
      Normal_Vert(:,contador_indices+13)=Nor_i
      Normal_Vert(:,contador_indices+14)=Nor3
      Tri(:,(count-1)*4+1)=contador_indices*[1, 1, 1, 1, 1, 1]+[0, 2, 9, 1, 6, 5]
      Tri(:,(count-1)*4+2)=contador_indices*[1, 1, 1, 1, 1, 1]+[2, 11, 9, 7, 10, 6]
      Tri(:,(count-1)*4+3)=contador_indices*[1, 1, 1, 1, 1, 1]+[2, 4, 11, 3, 8, 7]
      Tri(:,(count-1)*4+4)=contador_indices*[1, 1, 1, 1, 1, 1]+[9, 11, 14, 10, 13, 12]
      contador_indices=contador_indices+15
    enddo
    Geometry1%npoints=Geometry1%ntri*15
    Geometry1%ntri=Geometry1%ntri*4
    m=Geometry1%npoints
    N=Geometry1%ntri
    Geometry1%n_Sf_points=N*n_order_sf

    if (allocated(Geometry1%Points)) then
      deallocate(Geometry1%Points)
    endif
    allocate(Geometry1%Points(3,m))
    if (allocated(Geometry1%Tri)) then
      deallocate(Geometry1%Tri)
    endif
    allocate(Geometry1%Tri(6,N))
    if (allocated(Geometry1%Normal_Vert)) then
      deallocate(Geometry1%Normal_Vert)
    endif
    allocate(Geometry1%Normal_Vert(3,m))

    Geometry1%Points=Points
    Geometry1%Tri=Tri
    Geometry1%Normal_Vert=Normal_Vert

    deallocate(Points)
    deallocate(Normal_Vert)
    deallocate(Tri)

    return
  end subroutine refine_geometry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine funcion_Base_Points(Geometry1)
    implicit none

    !List of calling arguments
    type (Geometry), intent(inout) :: Geometry1

    !List of local variables
    double precision U(Geometry1%n_order_sf),V(Geometry1%n_order_sf),w(Geometry1%n_order_sf)
    double precision UV(2,Geometry1%n_order_sf)
    double precision P1(3),P2(3),P3(3),P4(3),P5(3),P6(3),N1(3),N2(3),N3(3)
    double precision F_x(Geometry1%n_order_sf),F_y(Geometry1%n_order_sf)
    double precision F_z(Geometry1%n_order_sf),dS(Geometry1%n_order_sf)
    double precision nP_x(Geometry1%n_order_sf),nP_y(Geometry1%n_order_sf),nP_z(Geometry1%n_order_sf)
    double precision U_x(Geometry1%n_order_sf),U_y(Geometry1%n_order_sf),U_z(Geometry1%n_order_sf)
    double precision V_x(Geometry1%n_order_sf),V_y(Geometry1%n_order_sf),V_z(Geometry1%n_order_sf)

    integer :: count,n_order_sf, itype, npols, norder_smooth,npol_smooth
    double precision, allocatable :: umatr(:,:),vmatr(:,:)


    
    norder_smooth = Geometry1%norder_smooth
    n_order_sf=Geometry1%n_order_sf

    npol_smooth = (norder_smooth+1)*(norder_smooth+2)/2
    allocate(umatr(npol_smooth,npol_smooth))
    allocate(vmatr(npol_smooth,npol_smooth))
    

    !
    ! get the smooth nodes and quadrature weights
    !
     npols = npol_smooth
     call vioreanu_simplex_quad(norder_smooth,npols,UV,umatr,vmatr,w)
     U = UV(1,:)
     V = UV(2,:)
    
    


    
    if (allocated(Geometry1%Base_Points)) then
      deallocate(Geometry1%Base_Points)
    endif
    if (allocated(Geometry1%Base_Points_N)) then
      deallocate(Geometry1%Base_Points_N)
    endif
    if (allocated(Geometry1%Base_Points_U)) then
      deallocate(Geometry1%Base_Points_U)
    endif
    if (allocated(Geometry1%Base_Points_V)) then
      deallocate(Geometry1%Base_Points_V)
    endif
    if (allocated(Geometry1%w_smooth)) then
      deallocate(Geometry1%w_smooth)
    endif
    allocate(Geometry1%Base_Points(3,Geometry1%n_Sf_points))
    allocate(Geometry1%Base_Points_N(3,Geometry1%n_Sf_points))
    allocate(Geometry1%Base_Points_U(3,Geometry1%n_Sf_points))
    allocate(Geometry1%Base_Points_V(3,Geometry1%n_Sf_points))
    allocate(Geometry1%w_smooth(Geometry1%n_Sf_points))
    do count=1,Geometry1%ntri
      P1=Geometry1%Points(:,Geometry1%Tri(1,count))
      P2=Geometry1%Points(:,Geometry1%Tri(2,count))
      P3=Geometry1%Points(:,Geometry1%Tri(3,count))
      P4=Geometry1%Points(:,Geometry1%Tri(4,count))
      P5=Geometry1%Points(:,Geometry1%Tri(5,count))
      P6=Geometry1%Points(:,Geometry1%Tri(6,count))
      N1=Geometry1%Normal_Vert(:,Geometry1%Tri(1,count))
      N2=Geometry1%Normal_Vert(:,Geometry1%Tri(2,count))
      N3=Geometry1%Normal_Vert(:,Geometry1%Tri(3,count))
      call eval_quadratic_patch_UV(P1,P2,P3,P4,P5,P6,U,V,F_x,F_y,F_z,U_x,U_y,U_z,V_x,V_y,V_z,nP_x,nP_y,nP_z,dS,n_order_sf)
      nP_x=N1(1)*(1.0d0-U-V)+N2(1)*U+N3(1)*V
      nP_y=N1(2)*(1.0d0-U-V)+N2(2)*U+N3(2)*V
      nP_z=N1(3)*(1.0d0-U-V)+N2(3)*U+N3(3)*V
      Geometry1%Base_Points(1,(count-1)*n_order_sf+1:(count)*n_order_sf)=F_x
      Geometry1%Base_Points(2,(count-1)*n_order_sf+1:(count)*n_order_sf)=F_y
      Geometry1%Base_Points(3,(count-1)*n_order_sf+1:(count)*n_order_sf)=F_z
      Geometry1%Base_Points_N(1,(count-1)*n_order_sf+1:(count)*n_order_sf)=nP_x
      Geometry1%Base_Points_N(2,(count-1)*n_order_sf+1:(count)*n_order_sf)=nP_y
      Geometry1%Base_Points_N(3,(count-1)*n_order_sf+1:(count)*n_order_sf)=nP_z
      Geometry1%Base_Points_U(1,(count-1)*n_order_sf+1:(count)*n_order_sf)=U_x
      Geometry1%Base_Points_U(2,(count-1)*n_order_sf+1:(count)*n_order_sf)=U_y
      Geometry1%Base_Points_U(3,(count-1)*n_order_sf+1:(count)*n_order_sf)=U_z

      Geometry1%Base_Points_V(1,(count-1)*n_order_sf+1:(count)*n_order_sf)=V_x
      Geometry1%Base_Points_V(2,(count-1)*n_order_sf+1:(count)*n_order_sf)=V_y
      Geometry1%Base_Points_V(3,(count-1)*n_order_sf+1:(count)*n_order_sf)=V_z

      Geometry1%w_smooth((count-1)*n_order_sf+1:(count)*n_order_sf)=w
    enddo
    return
  end subroutine funcion_Base_Points

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  
subroutine find_smooth_surface(Geometry1, Feval_stuff_1, adapt_flag, ier)
  implicit none

  !
  ! this is the main driver routine that finds the smooth surface via
  ! convolution
  !
  double precision :: my_cross(3)

  !List of calling arguments
  type (Geometry), intent(inout) :: Geometry1
  type ( Feval_stuff ), pointer :: Feval_stuff_1
  integer, intent(in) :: adapt_flag

  !List of local variables
  integer :: flag,maxiter,ipoint,count,n_order_sf,itri
  double precision :: h(Geometry1%n_Sf_points),tol
  !double precision :: Norm_N
  double precision :: r_t(3,Geometry1%n_Sf_points)
  double precision :: F(Geometry1%n_Sf_points)
  double precision :: grad_F(3,Geometry1%n_Sf_points)
  double precision :: dVdu(3),dVdv(3),dhdu,dhdv,h_var
  integer ier

  ier = 0

  !
  ! check the memory allocation
  !
  
  if (allocated(Geometry1%du_smooth)) then
    deallocate(Geometry1%du_smooth)
  endif
  allocate(Geometry1%du_smooth(3,Geometry1%n_Sf_points))
        
  if (allocated(Geometry1%dv_smooth)) then
    deallocate(Geometry1%dv_smooth)
  endif
  allocate(Geometry1%dv_smooth(3,Geometry1%n_Sf_points))
        
  if (allocated(Geometry1%ds_smooth)) then
    deallocate(Geometry1%ds_smooth)
  endif
  allocate(Geometry1%ds_smooth(Geometry1%n_Sf_points))

  if (allocated(Geometry1%S_smooth)) then
    deallocate(Geometry1%S_smooth)
  endif
  allocate(Geometry1%S_smooth(3,Geometry1%n_Sf_points))

  if (allocated(Geometry1%N_smooth)) then
    deallocate(Geometry1%N_smooth)
  endif
  allocate(Geometry1%N_smooth(3,Geometry1%n_Sf_points))

  if (allocated(Geometry1%ru_smooth)) then
    deallocate(Geometry1%ru_smooth)
  endif
  allocate(Geometry1%ru_smooth(3,Geometry1%n_Sf_points))

  if (allocated(Geometry1%rv_smooth)) then
    deallocate(Geometry1%rv_smooth)
  endif
  allocate(Geometry1%rv_smooth(3,Geometry1%n_Sf_points))



   if (.not.allocated(Geometry1%height)) then
     allocate (Geometry1%height(Geometry1%n_Sf_points))

     do count=1,Geometry1%n_Sf_points
       Geometry1%height(count)=0.0d0
       h(count)=0.0d0
     enddo
   else

     do count=1,Geometry1%n_Sf_points
       h(count)=Geometry1%height(count)
     enddo
   endif

  !
  ! set some parameters for the newton routine for finding the surface
  !
  tol = 1.0d-9
  maxiter = 14
  flag = 0


  call My_Newton(h, tol, maxiter, Geometry1, flag, Feval_stuff_1, &
      adapt_flag, grad_F, r_t)
  
  if (flag .eq. 1) then
    ier = 3
    return
  end if

  if (flag.eq. 3) then
     ier = 4
     return
  endif

  !
  ! copy over the points on the surface
  !
  do count = 1,Geometry1%n_Sf_points
    Geometry1%S_smooth(:,count) = r_t(:,count)
  enddo

  !
  ! compute various partial derivatives
  !
  n_order_sf = Geometry1%n_order_sf
  
  do count=1,Geometry1%n_Sf_points

    Geometry1%N_smooth(:,count) = -1.0d0*grad_F(:,count) &
        /(sqrt(grad_F(1,count)**2 + grad_F(2,count)**2 &
        + grad_F(3,count)**2))

    itri = (count-1)/n_order_sf+1
    h_var=h(count)

    dVdu = Geometry1%Normal_Vert(:,Geometry1%Tri(2,itri)) &
        - Geometry1%Normal_Vert(:,Geometry1%Tri(1,itri))

    dVdv = Geometry1%Normal_Vert(:,Geometry1%Tri(3,itri)) &
        - Geometry1%Normal_Vert(:,Geometry1%Tri(1,itri))

    dhdu = -dot_product(Geometry1%N_smooth(:,count), &
        dVdu*h_var+Geometry1%Base_Points_U(:,count)) &
        /dot_product(Geometry1%N_smooth(:,count) &
        ,Geometry1%Base_Points_N(:,count))

    dhdv = -dot_product(Geometry1%N_smooth(:,count), &
        dVdv*h_var+Geometry1%Base_Points_V(:,count)) &
        /dot_product(Geometry1%N_smooth(:,count), &
        Geometry1%Base_Points_N(:,count))

    Geometry1%ru_smooth(:,count) = dVdu*h_var &
        + Geometry1%Base_Points_N(:,count)*dhdu &
        + Geometry1%Base_Points_U(:,count)

    Geometry1%rv_smooth(:,count) = dVdv*h_var &
        + Geometry1%Base_Points_N(:,count)*dhdv &
        + Geometry1%Base_Points_V(:,count)

    Geometry1%du_smooth(:,count)=Geometry1%ru_smooth(:,count)
    Geometry1%dv_smooth(:,count)=Geometry1%rv_smooth(:,count)
    call crossproduct(Geometry1%ru_smooth(:,count)&
      &,Geometry1%rv_smooth(:,count),my_cross)
    Geometry1%ds_smooth(count)=norm2(my_cross)

    call crossproduct(Geometry1%ru_smooth(:,count), &
        Geometry1%rv_smooth(:,count), my_cross)

    Geometry1%w_smooth(count) = &
        Geometry1%w_smooth(count)*norm2(my_cross)

    Geometry1%ru_smooth(:,count) = &
        Geometry1%ru_smooth(:,count)/norm2(Geometry1%ru_smooth(:,count))

    call crossproduct(Geometry1%N_smooth(:,count), &
        Geometry1%ru_smooth(:,count), my_cross)

    Geometry1%rv_smooth(:,count) = my_cross

  enddo

  !
  ! copy over the h height function
  !
  do count=1,Geometry1%n_Sf_points
    Geometry1%height(count)=h(count)
  enddo

  return
end subroutine find_smooth_surface






subroutine My_Newton(x,tol,maxiter,Geometry1,flag, &
    Feval_stuff_1 ,adapt_flag, grad_F, r_t)
  implicit none

  !List of calling arguments
  type (Geometry), intent(in) :: Geometry1
  integer, intent(in) :: maxiter
  double precision, intent(inout) :: x(Geometry1%n_Sf_points)
  double precision, intent(in) :: tol
  integer, intent(out) :: flag
  type ( Feval_stuff ), pointer :: Feval_stuff_1      !! data type that
  integer, intent(in) :: adapt_flag
  double precision, intent(inout) :: r_t(3,Geometry1%n_Sf_points)
  double precision, intent(inout) :: grad_F(3,Geometry1%n_Sf_points)


  !
  ! this routine runs a newton iteration to find the smooth surface
  !
  ! Input:
  !   x - initial guess at the h function
  !   tol - pointwise convergence criterion for Newton
  !   maxiter - max iteration count
  !   Geometry1 - all the geo info
  !   Feval_stuff_1 -
  !   adapt_flag - determines what kind of smoothing kernel to use
  !
  ! Output:
  !   x - the h function from the surface
  !   flag - whether or not Newton succeeded
  !   grad_F - gradient of SOMETHING
  !   r_t - points on the smooth surface
  !   
  
  !List of local variables
  integer count,count2
  double precision F(Geometry1%n_Sf_points)
  double precision :: dF(Geometry1%n_Sf_points), err(Geometry1%n_Sf_points)
  integer  flag_con(Geometry1%n_Sf_points)

  ! initialize the errors and convergence flags
  do count2 = 1,Geometry1%n_Sf_points
    err(count2) = tol+1.0d0
    flag_con(count2) = 0
  enddo

  count=1
  flag=0

 4999 format(i10,i9,5x,e10.2)

  ! print out convergence information
  write (*,*) 'iteration,  #targets,   err'  
  write (*,4999) count, &
        Geometry1%n_Sf_points-sum(flag_con), 1.0d0 

  ! the x function has been initialized to 0
  do while ( (maxval(err)>tol) .and. (count<maxiter) )

    ! call a step of Newton

    !print *, "Entering fun_roots_derivative"

    call fun_roots_derivative(x, Geometry1, F, dF, Feval_stuff_1, &
        adapt_flag, flag_con, grad_F, r_t)

    count = count+1
    do count2 = 1,Geometry1%n_Sf_points

      if (flag_con(count2) .eq. 0) then

        err(count2) = F(count2)/dF(count2)
        x(count2) = x(count2) - err(count2)!    Mucho ojo con esto

        err(count2) = abs(err(count2))
        
        if ( abs(x(count2)) .lt. 1.0d0) then
          err(count2) = abs(err(count2))
        else
          err(count2) = abs(err(count2)/x(count2))
        end if
        
        
        if (abs(err(count2))<tol) then
          flag_con(count2)=1
        endif
      endif
    enddo

    write (*,4999) count, &
        Geometry1%n_Sf_points-sum(flag_con), maxval(err)
    if (maxval(err) > 1.1d0) then
       flag = 3
       return
    endif
        
  end do

!  stop

  
  if ( maxval(err)>tol ) then
    flag=1
    return
  endif

  return
end subroutine My_Newton






  
  subroutine check_Gauss(Geometry1,x0,y0,z0,err_rel)
    implicit none

    !List of calling arguments
    type (Geometry), intent(inout) :: Geometry1
    double precision, intent(in) :: x0,y0,z0
    double precision, intent(out) :: err_rel

    !List of local variables
    integer umio,count1,count2,flag,n_order_sf
    double precision  F,Ex,Ey,Ez,R,x,y,z,pi,w,nx,ny,nz

    pi=3.141592653589793238462643383d0
    F=0.0d0

    do count1=1,Geometry1%n_Sf_points
      x=Geometry1%S_smooth(1,count1)
      y=Geometry1%S_smooth(2,count1)
      z=Geometry1%S_smooth(3,count1)
      w=Geometry1%w_smooth(count1)
      nx=Geometry1%N_smooth(1,count1)
      ny=Geometry1%N_smooth(2,count1)
      nz=Geometry1%N_smooth(3,count1)

      R=sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
      Ex=(x-x0)/(4*pi*R**3)
      Ey=(y-y0)/(4*pi*R**3)
      Ez=(z-z0)/(4*pi*R**3)
      F=F+(Ex*nx+Ey*ny+Ez*nz)*w
    enddo
    err_rel=abs(F-1)
    
    return
  end subroutine check_Gauss





  
  subroutine fun_roots_derivative(h, Geometry1, F, dF, &
      Feval_stuff_1, &
      adapt_flag, flag_con, grad_F, r_t)
    implicit none

    !List of calling arguments
    type (Geometry), intent(in) :: Geometry1
    double precision, intent(in) :: h(Geometry1%n_Sf_points)
    double precision, intent(out) ::  F(Geometry1%n_Sf_points), &
        dF(Geometry1%n_Sf_points)
    type ( Feval_stuff ), pointer :: Feval_stuff_1
    
    integer, intent(in) :: adapt_flag
    integer, intent(in) :: flag_con(Geometry1%n_Sf_points)
    double precision, intent(inout) :: r_t(3,Geometry1%n_Sf_points), &
        grad_F(3,Geometry1%n_Sf_points)



    !List of local variables
    !double precision Norm_N
    integer :: count, ipointer, ntarg
    double precision, allocatable :: r_t2(:,:),v_norm(:,:),grad_F2(:,:),F2(:)

    allocate(r_t2(3,Geometry1%n_Sf_points-sum(flag_con)))
    allocate(v_norm(3,Geometry1%n_Sf_points-sum(flag_con)))
    allocate(grad_F2(3,Geometry1%n_Sf_points-sum(flag_con)))
    allocate(F2(Geometry1%n_Sf_points-sum(flag_con)))


    ipointer=1
    do count=1,Geometry1%n_Sf_points

      ! compute the current estimate for the point on the surface
      r_t(:,count) = Geometry1%Base_Points(:,count) &
          + Geometry1%Base_Points_N(:,count)*h(count)

      ! if not converged, add it to the list of targets
      if (flag_con(count) .eq. 0) then
        r_t2(:,ipointer) = r_t(:,count)
        v_norm(:,ipointer) = Geometry1%Base_Points_N(:,count)
        ipointer = ipointer+1
      endif
    enddo

    !
    ! get values and gradients
    !
    ntarg = Geometry1%n_Sf_points-sum(flag_con)

    call eval_density_grad_FMM(Geometry1, r_t2, v_norm, &
        ntarg, F2, grad_F2, Feval_stuff_1, adapt_flag)

    !
    ! update the points
    !
    ipointer=1
    do count=1,Geometry1%n_Sf_points

      if ( flag_con(count) .eq. 0) then

        dF(count)=(grad_F2(1,ipointer)*Geometry1%Base_Points_N(1,count)+grad_F2(2,ipointer)*&
            &Geometry1%Base_Points_N(2,count)+grad_F2(3,ipointer)*Geometry1%Base_Points_N(3,count))
        F(count)=F2(ipointer)
        r_t(:,count) = r_t2(:,ipointer)
        grad_F(:,count) = grad_F2(:,ipointer)
        ipointer=ipointer+1
      endif

    enddo

    deallocate(r_t2)
    deallocate(v_norm)
    deallocate(grad_F2)
    deallocate(F2)

    return
  end subroutine fun_roots_derivative





  
  subroutine refine_geometry_smart(Geometry1)
    implicit none

    !List of calling arguments
    type (Geometry), intent(inout) :: Geometry1

    !List of local variables
    character ( len=100 ) plot_name
    integer count,contador_indices
    double precision, allocatable :: Points(:,:), Normal_Vert(:,:)
    integer, allocatable :: Tri(:,:)
    double precision, allocatable :: h_new(:),h_tri(:),h_1(:),h_2(:),h_3(:)
    double precision, allocatable :: h_4(:)

    double precision P1(3),P2(3),P3(3),P4(3),P5(3),P6(3)
    double precision Pa(3),Pb(3),Pc(3),Pd(3),Pe(3),Pf(3),Pg(3),Ph(3),Pi(3)
    double precision Nor1(3),Nor2(3),Nor3(3),Nor4(3),Nor5(3),Nor6(3)
    double precision Nor_a(3),Nor_b(3),Nor_c(3),Nor_d(3),Nor_e(3),Nor_f(3),Nor_g(3),Nor_h(3),Nor_i(3)
    double precision U(9),V(9)
    double precision F_x(9),F_y(9),F_z(9),dS(9)
    double precision nP_x(9),nP_y(9),nP_z(9)
    double precision U_x(9),U_y(9),U_z(9),V_x(9),V_y(9),V_z(9)
    double precision, allocatable :: ximat(:,:)
    integer m,N,n_order_aux,n_order_sf,istart
    integer nover,norder, itmp

    n_order_sf=Geometry1%n_order_sf

    allocate(h_tri(n_order_sf))
    allocate(h_1(n_order_sf))
    allocate(h_2(n_order_sf))
    allocate(h_3(n_order_sf))
    allocate(h_4(n_order_sf))

    nover = n_order_sf*4
    allocate(ximat(nover,n_order_sf))

    norder = Geometry1%norder_smooth
    call get_refine_interp_mat(norder,n_order_sf,nover,ximat)

    allocate(Points(3,Geometry1%ntri*15))
    allocate(Normal_Vert(3,Geometry1%ntri*15))
    allocate(Tri(6,Geometry1%ntri*4))
    allocate(h_new(Geometry1%n_Sf_points*4))

    contador_indices=1
    n_order_aux=9
    U=[0.2500d0,0.7500d0,0d0,0.2500d0,0.5000d0,0.7500d0,0.2500d0,0d0,0.2500d0]
    V=[0d0,0d0,0.2500d0,0.2500d0,0.2500d0,0.2500d0,0.5000d0,0.7500d0,0.7500d0]

    do count=1,Geometry1%ntri

      P1=Geometry1%Points(:,Geometry1%Tri(1,count))
      P2=Geometry1%Points(:,Geometry1%Tri(2,count))
      P3=Geometry1%Points(:,Geometry1%Tri(3,count))
      P4=Geometry1%Points(:,Geometry1%Tri(4,count))
      P5=Geometry1%Points(:,Geometry1%Tri(5,count))
      P6=Geometry1%Points(:,Geometry1%Tri(6,count))

      call eval_quadratic_patch_UV(P1,P2,P3,P4,P5,P6,U,V,F_x,F_y,F_z,U_x,U_y,U_z,V_x,V_y,V_z,nP_x,nP_y,nP_z,dS,n_order_aux)
      Pa=[F_x(1),F_y(1),F_z(1)]
      Pb=[F_x(2),F_y(2),F_z(2)]
      Pc=[F_x(3),F_y(3),F_z(3)]
      Pd=[F_x(4),F_y(4),F_z(4)]
      Pe=[F_x(5),F_y(5),F_z(5)]
      Pf=[F_x(6),F_y(6),F_z(6)]
      Pg=[F_x(7),F_y(7),F_z(7)]
      Ph=[F_x(8),F_y(8),F_z(8)]
      Pi=[F_x(9),F_y(9),F_z(9)]
      Nor1=Geometry1%Normal_Vert(:,Geometry1%Tri(1,count))
      Nor2=Geometry1%Normal_Vert(:,Geometry1%Tri(2,count))
      Nor3=Geometry1%Normal_Vert(:,Geometry1%Tri(3,count))
      Nor4=(Nor1+Nor2)/2.0d0
      Nor5=(Nor2+Nor3)/2.0d0
      Nor6=(Nor3+Nor1)/2.0d0
      Nor_a=(Nor1+Nor4)/2.0d0
      Nor_b=(Nor4+Nor2)/2.0d0
      Nor_c=(Nor1+Nor6)/2.0d0
      Nor_d=(Nor4+Nor6)/2.0d0
      Nor_e=(Nor4+Nor5)/2.0d0
      Nor_f=(Nor5+Nor2)/2.0d0
      Nor_g=(Nor5+Nor6)/2.0d0
      Nor_h=(Nor3+Nor6)/2.0d0
      Nor_i=(Nor3+Nor5)/2.0d0
      Points(:,contador_indices)=P1
      Points(:,contador_indices+1)=Pa
      Points(:,contador_indices+2)=P4
      Points(:,contador_indices+3)=Pb
      Points(:,contador_indices+4)=P2
      Points(:,contador_indices+5)=Pc
      Points(:,contador_indices+6)=Pd
      Points(:,contador_indices+7)=Pe
      Points(:,contador_indices+8)=Pf
      Points(:,contador_indices+9)=P6
      Points(:,contador_indices+10)=Pg
      Points(:,contador_indices+11)=P5
      Points(:,contador_indices+12)=Ph
      Points(:,contador_indices+13)=Pi
      Points(:,contador_indices+14)=P3

      Normal_Vert(:,contador_indices)=Nor1
      Normal_Vert(:,contador_indices+1)=Nor_a
      Normal_Vert(:,contador_indices+2)=Nor4
      Normal_Vert(:,contador_indices+3)=Nor_b
      Normal_Vert(:,contador_indices+4)=Nor2
      Normal_Vert(:,contador_indices+5)=Nor_c
      Normal_Vert(:,contador_indices+6)=Nor_d
      Normal_Vert(:,contador_indices+7)=Nor_e
      Normal_Vert(:,contador_indices+8)=Nor_f
      Normal_Vert(:,contador_indices+9)=Nor6
      Normal_Vert(:,contador_indices+10)=Nor_g
      Normal_Vert(:,contador_indices+11)=Nor5
      Normal_Vert(:,contador_indices+12)=Nor_h
      Normal_Vert(:,contador_indices+13)=Nor_i
      Normal_Vert(:,contador_indices+14)=Nor3
      Tri(:,(count-1)*4+1)=contador_indices*[1, 1, 1, 1, 1, 1]+[0, 2, 9, 1, 6, 5]
      Tri(:,(count-1)*4+2)=contador_indices*[1, 1, 1, 1, 1, 1]+[11, 9, 2, 10, 6, 7]![2, 11, 9, 7, 10, 6]
      Tri(:,(count-1)*4+3)=contador_indices*[1, 1, 1, 1, 1, 1]+[2, 4, 11, 3, 8, 7]
      Tri(:,(count-1)*4+4)=contador_indices*[1, 1, 1, 1, 1, 1]+[9, 11, 14, 10, 13, 12]

      h_tri=Geometry1%height((count-1)*n_order_sf+1:(count)*n_order_sf)

      istart = (count-1)*nover+1
      call dmatvec(nover,n_order_sf,ximat,h_tri,h_new(istart))


      contador_indices=contador_indices+15

    enddo

    Geometry1%npoints=Geometry1%ntri*15
    Geometry1%ntri=Geometry1%ntri*4
    m=Geometry1%npoints
    N=Geometry1%ntri
    Geometry1%n_Sf_points=N*n_order_sf

    if (allocated(Geometry1%Points)) then
      deallocate(Geometry1%Points)
    endif
    allocate(Geometry1%Points(3,m))

    if (allocated(Geometry1%Tri)) then
      deallocate(Geometry1%Tri)
    endif
    allocate(Geometry1%Tri(6,N))

    if (allocated(Geometry1%Normal_Vert)) then
      deallocate(Geometry1%Normal_Vert)
    endif
    allocate(Geometry1%Normal_Vert(3,m))

    Geometry1%Points = Points
    Geometry1%Tri = Tri
    Geometry1%Normal_Vert = Normal_Vert

    deallocate(Points)
    deallocate(Normal_Vert)
    deallocate(Tri)


    deallocate(Geometry1%height)

    allocate(Geometry1%height(Geometry1%n_Sf_points))
    Geometry1%height=h_new

    return
  end subroutine refine_geometry_smart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module Mod_Smooth_Surface

