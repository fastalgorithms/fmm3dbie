!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!                    MODULE MAIN INFORMATION                          !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! DESCRIPTION:
!!!
!!! This module contains functions and data types to allow the
!!! evaluation of the function F Notice that we have already
!!! substracted 0.5d0, therefore Newton has to find the zero set
!!! F(x,y,z) = 0.0d0
!!!
!!! INDEX OF SUBROUTINES AND FUNCTIONS:
!!!
!!! Public:
!!!
!!!  ! generate all the data structure required to evaluate F
!!!  ! including the trees for the FMM and sgma generator
!!!  1-start_Feval                        
!!!
!!!  ! function that evaluates F and gradient at arbitrary points in space.
!!!  2-eval_density_grad_FMM                
!!!
!!! Private: (none)
!!!
!!! INDEX OF DERIVED TYPES:
!!!
!!! Public:
!!!
!!!   Feval_stuff - All the information required to evaluate F,
!!!                  including Manas' trees, and the the tree required
!!!                  to evaluate sigma
!!!
!!! Private: (none)
!!!
!!!
!!! MODULES REQUIRED:
!!!
!!! 1-Mod_Fast_Sigma
!!! 2-ModType_Smooth_Surface
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module Mod_Feval

  use Mod_TreeLRD
  use Mod_Fast_Sigma
  use ModType_Smooth_Surface

  implicit none
  
  !
  ! type definitions
  !
  type Feval_stuff
     type (  Fast_Sigma_stuff ), pointer :: FSS_1 => null ()
     type ( TreeLRD ), pointer :: Tree_local => null ()
     !!All the stuff for Manas' tree
     integer ltree,norder,nlevels,nboxes,nt2
     integer, allocatable :: itree(:)
     integer, allocatable :: iptr(:)
     real ( kind = 8 ) eps
     real ( kind = 8 ), allocatable :: treecenters(:,:)
     real ( kind = 8 ), allocatable :: boxsize(:)
     real ( kind = 8 ), allocatable :: fcoeffs(:)
     real ( kind = 8 ), allocatable :: fcoeffsx(:)
     real ( kind = 8 ), allocatable :: fcoeffsy(:)
     real ( kind = 8 ), allocatable :: fcoeffsz(:)

  end type Feval_stuff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!    END DEFINITION OF TYPES     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!    DEFINITION OF SUBROUTINES     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !public :: start_Feval
  !public :: eval_density_grad_FMM

contains




  



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




  subroutine start_Feval_tree(Feval_stuff_1,Geometry1,rlam)
    implicit none

    !
    ! This function initializes all the data structure required to
    ! evaluate the F function (including FMM trees and sgma trees)
    !

    !List of calling arguments
    type ( Feval_stuff ), pointer :: Feval_stuff_1      !! data type that contains all the information to evaluate F
    type ( Geometry ), intent(inout)  :: Geometry1      !! data type that contains all the information about the geometry

    !List of local variables
    integer count1,count,n_targets

    real *8 rlam

    ! Allocate all arrays required
    if (.not.(associated(Feval_stuff_1))) then
      allocate(Feval_stuff_1)
      call setup_tree_sigma_geometry(Feval_stuff_1%FSS_1,Geometry1,rlam)   !! Setup the tree structure to evaluate sgma function
    endif

    return
  end subroutine start_Feval_tree
  
  



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine eval_density_grad_FMM(Geometry1, targets, v_norm, &
      n_targets, F,grad_F,Fev_stf_1,adapt_flag)
    implicit none
    !
    ! This function evaluates the F function (and gradient) using the FMM and the sgma evaluator
    !! Flags to select the different versions of F
    !integer, intent(in) :: adapt_flag,n_targets
    integer  :: adapt_flag
    integer  :: n_targets

    !List of calling arguments
    type ( Geometry ), intent(in) :: Geometry1
    real ( kind = 8 ) , intent(in) :: targets(3,n_targets)
    real ( kind = 8 ) , intent(in) :: v_norm(3,n_targets)
    real ( kind = 8 ), intent(out) ::  F(n_targets), grad_F(3,n_targets)

    !! all the information required to evaluate F (trees..)
    type ( Feval_stuff ), pointer :: Fev_stf_1    


    !List of local variables
    !! Different flags to run appropriately the FMM
    integer ier, iprec,ifcharge,ifdipole,ifpottarg,iffldtarg   
    integer count,count2,count1
    integer n_sources
    real ( kind = 8 ) tfmm
    real ( kind = 8 ), allocatable :: sgma(:),sgma_grad(:,:)
    integer , allocatable :: flag_error(:)
    real ( kind = 8 ), allocatable :: missed_Points(:,:)
    integer ipointer
    character (len=100) plot_name
    double precision :: t0, t1, telap, sgma_max, eps
    
    !$ double precision :: omp_get_wtime

    allocate(sgma(n_targets))
    allocate(sgma_grad(3,n_targets))
    !
    ! obtain the value of sgma at the target points, needed to make
    ! the FMM call
    !
    !write (*,*) 'START eval sigma',n_targets

    call cpu_time(t0)
!$   t0 = omp_get_wtime()    
    call function_eval_sigma(Fev_stf_1%FSS_1, targets, n_targets, sgma, &
        sgma_grad(1,:), sgma_grad(2,:),sgma_grad(3,:), adapt_flag)
    call cpu_time(t1)
!$   t1 = omp_get_wtime()
    telap = t1-t0
    !print *, 'time for eval sigma = ', telap

    !       read (*,*)
    !write (*,*) 'STOP eval sigma'
    


    ifcharge=0
    ifdipole=1
    ifpottarg=1
    iffldtarg=1
    iprec=3
    ier=0
    n_sources=Geometry1%n_Sk_points

    !! FMM call (this includes the target dependent local corrections with the erf function)
    !    write (*,*) 'target: ',targets(:,1)


    call cpu_time(t0)
    !$ t0 = omp_get_wtime()
    
    eps = 0.51d-9

    call tfmm3d_setsub(eps, n_sources, Geometry1%skeleton_Points, &
      Geometry1%skeleton_N, Geometry1%skeleton_w, n_targets, &
      targets, sgma, sgma_grad, F, grad_F)

    call cpu_time(t1)

    !$ t1 = omp_get_wtime()
    telap = t1-t0


        !
        ! compute the error in the levelset function
        !
    do count2=1,n_targets
      F(count2)=F(count2)-0.5d0
    enddo



    deallocate(sgma)
    deallocate(sgma_grad)


    return
  end subroutine eval_density_grad_FMM




  
  subroutine generate_dummy_targets(Geometry1,Fev_stf_1,adapt_flag)
    implicit none
    ! This is to generate the dummy targets

    !List of calling arguments
    type ( Geometry ), intent(inout) :: Geometry1                   !! data type that contains all the information about the geometry
    type ( Feval_stuff ), pointer :: Fev_stf_1      !! data type that contains all the information to evaluate F
    integer, intent(in) :: adapt_flag

    !List of local variables
    real ( kind = 8 ), allocatable :: dummy_t_1(:,:),dummy_t_2(:,:)
    real ( kind = 8 ), allocatable :: sgma(:), sgma_grad(:,:)
    integer count1,count2,ipcount,n_order,n_tri_dummy,n_ref
    character (len=100) plot_name
    real ( kind = 8 ) hhh


    n_order=Geometry1%n_Sf_points/Geometry1%ntri
    allocate(dummy_t_1(3,Geometry1%n_Sf_points))
    dummy_t_1=Geometry1%S_smooth
    allocate(dummy_t_2(3,Geometry1%n_Sf_points*4))
    n_ref=1 !!! MUCHO OJO CON ESTO; HAY QUE PONERLO A $ O ASÍ
    Geometry1%n_dummy_targ=Geometry1%n_Sf_points*(1-4**(n_ref+1))/(1-4)-Geometry1%n_Sf_points
    write (*,*) 'initial number of Dummy targets: ',Geometry1%n_dummy_targ
    allocate(Geometry1%Dummy_targ(3,Geometry1%n_dummy_targ))

    allocate(sgma(Geometry1%n_dummy_targ))
    allocate(sgma_grad(3,Geometry1%n_dummy_targ))

    ipcount=1
    do count1=1,n_ref
      call refine_mesh_1(dummy_t_1,n_order,Geometry1%ntri*4**(count1-1),dummy_t_2)
      deallocate(dummy_t_1)
      allocate(dummy_t_1(3,Geometry1%ntri*4**(count1)*n_order))
      dummy_t_1=dummy_t_2
      deallocate(dummy_t_2)
      if (count1<n_ref) then
        allocate(dummy_t_2(3,Geometry1%ntri*4**(count1+1)*n_order))
      endif
      Geometry1%Dummy_targ(:,ipcount:ipcount+Geometry1%ntri*4**(count1-1)*n_order-1)=dummy_t_1
      ipcount=ipcount+Geometry1%ntri*4**(count1-1)*n_order
    enddo


    call function_eval_sigma(Fev_stf_1%FSS_1,Geometry1%Dummy_targ,Geometry1%n_dummy_targ,sgma,&
        &sgma_grad(1,:),sgma_grad(2,:),sgma_grad(3,:),adapt_flag)

    if (allocated(dummy_t_2)) then
      deallocate(dummy_t_2)
    endif
    !allocate(dummy_t_2(3,Geometry1%n_dummy_targ*7))
    allocate(dummy_t_2(3,Geometry1%n_dummy_targ*8))


    hhh=1.0d-4
    do count1=1,Geometry1%n_dummy_targ
      dummy_t_2(1,count1)=Geometry1%Dummy_targ(1,count1)+sgma(count1)*hhh
      dummy_t_2(2,count1)=Geometry1%Dummy_targ(2,count1)+sgma(count1)*hhh
      dummy_t_2(3,count1)=Geometry1%Dummy_targ(3,count1)+sgma(count1)*hhh

      dummy_t_2(1,Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(1,count1)+sgma(count1)*hhh
      dummy_t_2(2,Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(2,count1)+sgma(count1)*hhh
      dummy_t_2(3,Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(3,count1)-sgma(count1)*hhh

      dummy_t_2(1,2*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(1,count1)+sgma(count1)*hhh
      dummy_t_2(2,2*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(2,count1)-sgma(count1)*hhh
      dummy_t_2(3,2*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(3,count1)+sgma(count1)*hhh

      dummy_t_2(1,3*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(1,count1)+sgma(count1)*hhh
      dummy_t_2(2,3*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(2,count1)-sgma(count1)*hhh
      dummy_t_2(3,3*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(3,count1)-sgma(count1)*hhh

      dummy_t_2(1,4*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(1,count1)-sgma(count1)*hhh
      dummy_t_2(2,4*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(2,count1)+sgma(count1)*hhh
      dummy_t_2(3,4*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(3,count1)+sgma(count1)*hhh

      dummy_t_2(1,5*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(1,count1)-sgma(count1)*hhh
      dummy_t_2(2,5*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(2,count1)+sgma(count1)*hhh
      dummy_t_2(3,5*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(3,count1)-sgma(count1)*hhh

      dummy_t_2(1,6*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(1,count1)-sgma(count1)*hhh
      dummy_t_2(2,6*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(2,count1)-sgma(count1)*hhh
      dummy_t_2(3,6*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(3,count1)+sgma(count1)*hhh


      dummy_t_2(1,7*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(1,count1)-sgma(count1)*hhh
      dummy_t_2(2,7*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(2,count1)-sgma(count1)*hhh
      dummy_t_2(3,7*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(3,count1)-sgma(count1)*hhh



      !dummy_t_2(1,count1)=Geometry1%Dummy_targ(1,count1)+sgma(count1)*hhh
      !dummy_t_2(2,count1)=Geometry1%Dummy_targ(2,count1)
      !dummy_t_2(3,count1)=Geometry1%Dummy_targ(3,count1)

      !dummy_t_2(1,Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(1,count1)-sgma(count1)*hhh
      !dummy_t_2(2,Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(2,count1)
      !dummy_t_2(3,Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(3,count1)

      !dummy_t_2(1,2*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(1,count1)
      !dummy_t_2(2,2*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(2,count1)+sgma(count1)*hhh
      !dummy_t_2(3,2*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(3,count1)

      !dummy_t_2(1,3*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(1,count1)
      !dummy_t_2(2,3*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(2,count1)-sgma(count1)*hhh
      !dummy_t_2(3,3*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(3,count1)

      !dummy_t_2(1,4*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(1,count1)
      !dummy_t_2(2,4*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(2,count1)
      !dummy_t_2(3,4*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(3,count1)+sgma(count1)*hhh

      !dummy_t_2(1,5*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(1,count1)
      !dummy_t_2(2,5*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(2,count1)
      !dummy_t_2(3,5*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(3,count1)-sgma(count1)*hhh

      !dummy_t_2(1,6*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(1,count1)
      !dummy_t_2(2,6*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(2,count1)
      !dummy_t_2(3,6*Geometry1%n_dummy_targ+count1)=Geometry1%Dummy_targ(3,count1)



    enddo

    deallocate(Geometry1%Dummy_targ)
    allocate(Geometry1%Dummy_targ(3,Geometry1%n_dummy_targ*7))
    Geometry1%Dummy_targ=dummy_t_2
    !        Geometry1%n_dummy_targ=Geometry1%n_dummy_targ*7
    Geometry1%n_dummy_targ=Geometry1%n_dummy_targ*8


    if (allocated(dummy_t_1)) then
      deallocate(dummy_t_1)
    endif
    if (allocated(dummy_t_2)) then
      deallocate(dummy_t_2)
    endif

    !        plot_name='./plot_tools/plot_dummy'
    !        call plot_curve_3D(Geometry1%Dummy_targ(1,:),Geometry1%Dummy_targ(2,:),Geometry1%Dummy_targ(3,:)&
    !         &,Geometry1%n_dummy_targ,plot_name)

    !   stop

    return
  end subroutine generate_dummy_targets

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine refine_mesh_1(Points_in,n_order,n_tri,Points_out)
    implicit none

    !List of calling arguments
    integer, intent(in) :: n_order,n_tri
    real ( kind = 8 ), intent(in) :: Points_in(3,n_tri*n_order)
    real ( kind = 8 ), intent(out) :: Points_out(3,n_tri*n_order*4)

    !List of local variables
    integer count,count1,count2,ipcount
    !    real ( kind = 8 ) sgma_tri(n_order)
    real ( kind = 8 ) x_tri(n_order),y_tri(n_order),z_tri(n_order)
    real ( kind = 8 ) x_tri_1(n_order),y_tri_1(n_order),z_tri_1(n_order)
    real ( kind = 8 ) x_tri_2(n_order),y_tri_2(n_order),z_tri_2(n_order)
    real ( kind = 8 ) x_tri_3(n_order),y_tri_3(n_order),z_tri_3(n_order)
    real ( kind = 8 ) x_tri_4(n_order),y_tri_4(n_order),z_tri_4(n_order)

    do count=1,n_tri
      x_tri=Points_in(1,(count-1)*n_order+1:(count)*n_order)
      y_tri=Points_in(2,(count-1)*n_order+1:(count)*n_order)
      z_tri=Points_in(3,(count-1)*n_order+1:(count)*n_order)
      if (n_order.eq.45) then
        call refine_tri45(x_tri,x_tri_1,x_tri_2,x_tri_3,x_tri_4)
        call refine_tri45(y_tri,y_tri_1,y_tri_2,y_tri_3,y_tri_4)
        call refine_tri45(z_tri,z_tri_1,z_tri_2,z_tri_3,z_tri_4)
      elseif (n_order.eq.78) then
        call refine_tri78(x_tri,x_tri_1,x_tri_2,x_tri_3,x_tri_4)
        call refine_tri78(y_tri,y_tri_1,y_tri_2,y_tri_3,y_tri_4)
        call refine_tri78(z_tri,z_tri_1,z_tri_2,z_tri_3,z_tri_4)
      else
        write (*,*) 'ERROR, n_order not found (different from 45 or 78); n_order = ', n_order
        read (*,*)
      endif
      Points_out(1,(count-1)*n_order*4+1:(count-1)*n_order*4+n_order)=x_tri_1
      Points_out(1,(count-1)*n_order*4+n_order+1:(count-1)*n_order*4+2*n_order)=x_tri_2
      Points_out(1,(count-1)*n_order*4+2*n_order+1:(count-1)*n_order*4+3*n_order)=x_tri_3
      Points_out(1,(count-1)*n_order*4+3*n_order+1:(count-1)*n_order*4+4*n_order)=x_tri_4

      Points_out(2,(count-1)*n_order*4+1:(count-1)*n_order*4+n_order)=y_tri_1
      Points_out(2,(count-1)*n_order*4+n_order+1:(count-1)*n_order*4+2*n_order)=y_tri_2
      Points_out(2,(count-1)*n_order*4+2*n_order+1:(count-1)*n_order*4+3*n_order)=y_tri_3
      Points_out(2,(count-1)*n_order*4+3*n_order+1:(count-1)*n_order*4+4*n_order)=y_tri_4

      Points_out(3,(count-1)*n_order*4+1:(count-1)*n_order*4+n_order)=z_tri_1
      Points_out(3,(count-1)*n_order*4+n_order+1:(count-1)*n_order*4+2*n_order)=z_tri_2
      Points_out(3,(count-1)*n_order*4+2*n_order+1:(count-1)*n_order*4+3*n_order)=z_tri_3
      Points_out(3,(count-1)*n_order*4+3*n_order+1:(count-1)*n_order*4+4*n_order)=z_tri_4
    enddo

    return
  end subroutine refine_mesh_1


  subroutine setup_edges_boundary(Geometry1)
    implicit none


    !List of calling arguments
    type ( Geometry ), intent(inout)  :: Geometry1

    !List of local variables
    integer, allocatable :: Edges_aux(:,:),Boundary_aux(:,:)
    integer icount, count1,count2,n_nodes_l
    integer, allocatable :: flags(:)
    real ( kind = 8 ) coef_x(3),coef_y(3),coef_z(3),P1(3),P2(3),P3(3),t(32),w(32)

    allocate(flags(Geometry1%npoints))
    allocate(Edges_aux(3,Geometry1%npoints))
    allocate(Geometry1%Boundary(3,Geometry1%ntri_sk))
    do count1=1,Geometry1%npoints
      flags(count1)=0
    enddo
    icount=1
    !    write (*,*) 'Empezamos'
    do count1=1,Geometry1%ntri_sk
      !        write (*,*) count1
      if (flags(Geometry1%Tri(4,count1)).eq.0) then
        !            write (*,*) '1 no existe'
        Edges_aux(1,icount)=Geometry1%Tri(1,count1)
        Edges_aux(2,icount)=Geometry1%Tri(4,count1)
        Edges_aux(3,icount)=Geometry1%Tri(2,count1)
        Geometry1%Boundary(1,count1)=icount
        flags(Geometry1%Tri(4,count1))=icount
        icount=icount+1
      else
        !            write (*,*) '1 ya existe'
        Geometry1%Boundary(1,count1)=-flags(Geometry1%Tri(4,count1))
      endif

      if (flags(Geometry1%Tri(5,count1)).eq.0) then
        Edges_aux(1,icount)=Geometry1%Tri(2,count1)
        Edges_aux(2,icount)=Geometry1%Tri(5,count1)
        Edges_aux(3,icount)=Geometry1%Tri(3,count1)
        Geometry1%Boundary(2,count1)=icount
        flags(Geometry1%Tri(5,count1))=icount
        icount=icount+1
      else
        Geometry1%Boundary(2,count1)=-flags(Geometry1%Tri(5,count1))
      endif

      if (flags(Geometry1%Tri(6,count1)).eq.0) then
        Edges_aux(1,icount)=Geometry1%Tri(3,count1)
        Edges_aux(2,icount)=Geometry1%Tri(6,count1)
        Edges_aux(3,icount)=Geometry1%Tri(1,count1)
        Geometry1%Boundary(3,count1)=icount
        flags(Geometry1%Tri(6,count1))=icount
        icount=icount+1
      else
        Geometry1%Boundary(3,count1)=-flags(Geometry1%Tri(6,count1))
      endif
    enddo
    icount=icount-1
    allocate(Geometry1%Edges(3,icount))
    Geometry1%Edges=Edges_aux(:,1:icount)
    n_nodes_l=32
    call Gauss1D(t,w,n_nodes_l)
    allocate(Geometry1%skelet_line_p(3,n_nodes_l,icount))
    allocate(Geometry1%skelet_line_dl(3,n_nodes_l,icount))
    do count1=1,icount
      P1=Geometry1%Points(:,Geometry1%Edges(1,count1))
      P2=Geometry1%Points(:,Geometry1%Edges(2,count1))
      P3=Geometry1%Points(:,Geometry1%Edges(3,count1))
      coef_x(1)=P1(1)
      coef_x(2)=-3*P1(1)+4*P2(1)-P3(1)
      coef_x(3)=2*P1(1)-4*P2(1)+2*P3(1)
      coef_y(1)=P1(2)
      coef_y(2)=-3*P1(2)+4*P2(2)-P3(2)
      coef_y(3)=2*P1(2)-4*P2(2)+2*P3(2)
      coef_z(1)=P1(3)
      coef_z(2)=-3*P1(3)+4*P2(3)-P3(3)
      coef_z(3)=2*P1(3)-4*P2(3)+2*P3(3)
      do count2=1,n_nodes_l
        Geometry1%skelet_line_p(1,count2,count1)=coef_x(1)+coef_x(2)*t(count2)+coef_x(3)*t(count2)**2
        Geometry1%skelet_line_p(2,count2,count1)=coef_y(1)+coef_y(2)*t(count2)+coef_y(3)*t(count2)**2
        Geometry1%skelet_line_p(3,count2,count1)=coef_z(1)+coef_z(2)*t(count2)+coef_z(3)*t(count2)**2
        Geometry1%skelet_line_dl(1,count2,count1)=(coef_x(2)+2.0d0*coef_x(3)*t(count2))*w(count2)
        Geometry1%skelet_line_dl(2,count2,count1)=(coef_y(2)+2.0d0*coef_y(3)*t(count2))*w(count2)
        Geometry1%skelet_line_dl(3,count2,count1)=(coef_z(2)+2.0d0*coef_z(3)*t(count2))*w(count2)
      enddo
    enddo
    !write (*,*) size(Geometry1%Points)
    !        read (*,*)
    return
  end subroutine setup_edges_boundary



  subroutine tpot_contour(source,dl,targ,v_norm,pot,fld)
    implicit none

    !List of calling arguments
    real ( kind = 8 ), intent(in) :: source(3),dl(3),targ(3),v_norm(3)
    real ( kind = 8 ), intent(out) ::  pot, fld(3)

    !List of local variables
    real ( kind = 8 ) F_aux(3)


    !    call Field_contour(source, v_norm, F_aux)
    call Field_contour(source-targ, v_norm, F_aux)

    pot=F_aux(1)*dl(1)+F_aux(2)*dl(2)+F_aux(3)*dl(3)
    !    write (*,*) dl


    return
  end subroutine tpot_contour

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tpot_contour2(source,dl,targ,n_vect_in,pot,fld)
    implicit none

    !List of calling arguments
    real ( kind = 8 ), intent(in) :: source(3),dl(3),targ(3),n_vect_in(3)
    real ( kind = 8 ), intent(out) ::  pot, fld(3)

    !List of local variables
    real ( kind = 8 ) F_aux(3),n_vect(3),norm_n,pt(3),f,dlxn(3),nxpt(3),pi
    real ( kind = 8 ) norm_pt,ptdotdlxn,v_aux1(3),v_aux2(3),norm_nxpt2

    pi=3.1415926535897932384626433d0
    norm_n=sqrt(n_vect_in(1)**2+n_vect_in(2)**2+n_vect_in(3)**2)
    n_vect=-1.0d0*n_vect_in/norm_n
    pt=source-targ
    norm_pt=sqrt(pt(1)**2+pt(2)**2+pt(3)**2)
    dlxn=cross_prod(dl,n_vect)
    nxpt=cross_prod(n_vect,pt)
    norm_nxpt2=nxpt(1)**2+nxpt(2)**2+nxpt(3)**2
    f=dot_prod(n_vect,pt)/norm_pt
    ptdotdlxn=dot_prod(pt,dlxn)
    pot=ptdotdlxn/norm_nxpt2*(1-f)
    pot=pot/(4.0d0*pi)
    v_aux1=-1.0d0*(n_vect/norm_pt-dot_prod(n_vect,pt)*pt/norm_pt**3)/norm_nxpt2
    v_aux2=2.0d0*(1-f)/norm_nxpt2**2*(cross_prod(n_vect,nxpt))
    fld=dlxn*(1-f)/norm_nxpt2+ptdotdlxn*(v_aux1+v_aux2)
    fld=fld/(4.0d0*pi)

    return
  end subroutine tpot_contour2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine Field_contour(point, vector, F)
    implicit none

    !List of calling arguments
    real ( kind = 8 ), intent(in) :: point(3), vector(3)
    real ( kind = 8 ), intent(out) :: F(3)

    !List of local variables
    real ( kind = 8 ) v_aux(3),v_aux2(3),normv,normv2,factor,normp,pi
    integer  count1, count2


    pi=3.1415926535897932384626433d0
    normp=sqrt(point(1)**2+point(2)**2+point(3)**2)
    normv=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
    v_aux2(1)=vector(1)/normv
    v_aux2(2)=vector(2)/normv
    v_aux2(3)=vector(3)/normv
    v_aux=point-v_aux2*dot_prod(v_aux2,point)
    normv2=(v_aux(1)**2+v_aux(2)**2+v_aux(3)**2)
    F=cross_prod(v_aux2,v_aux)/normv2
    factor=dot_prod(v_aux2,point)/normp
    F=F*(1.0d0-factor)/(4.0d0*pi)

    !        F(1)=-point(2)/2.0d0
    !        F(2)=point(1)/2.0d0
    !        F(3)=0.0d0

    !        F(1)=0.0d0
    !        F(2)=0.0d0
    !        F(3)=1.0d0

    return
  end subroutine Field_contour


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function cross_prod(v1,v2) result(v3)
    implicit none

    !List of calling arguments
    real ( kind = 8 ) v3(3)
    real ( kind = 8 ), intent(in) :: v1(3),v2(3)

    v3(1)=v1(2)*v2(3)-v1(3)*v2(2)
    v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
    v3(3)=v1(1)*v2(2)-v1(2)*v2(1)

  end function cross_prod


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function dot_prod(v1,v2) result(v3)
    implicit none

    !List of calling arguments
    real ( kind = 8 ) v3
    real ( kind = 8 ), intent(in) :: v1(3),v2(3)

    v3=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)

  end function dot_prod


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module Mod_Feval
