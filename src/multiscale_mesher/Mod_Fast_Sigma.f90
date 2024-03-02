!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!                    MODULE MAIN INFORMATION                          !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! DESCRIPTION:
!!!
!!! This module contains the sigma evaluator
!!! there are different versions that can be called with different flags
!!! there are different speeds (fast/brute force) that can be called with different flags
!!!
!!! INDEX OF SUBROUTINES AND FUNCTIONS:
!!!
!!! Public:
!!!
!!!  1-setup_tree_sigma_geometry          !! Initialize the tree, requires Geometry1
!!!  2-function_eval_sigma                !! Evaluate sigma and gradient of sigma at any point
!!!
!!! Private:
!!!
!!!  3-find_centroids_sigmas              !! Initialize centroids of the skeleton's triangles and
!!!                                       !! values of sigma at those points
!!!  4-function_eval_sigma_v0             !! computes sgma with no adaptivity at all (constant sgma
!!!                                       !! and derivatives=0)
!!!  5-fast_gaussian_global
!!!  6-fast_gaussian_box
!!!  7-fast_gaussian_box_grad
!!!  8-initial_guess_sgma                 !! provides a good initial guess for the value of sigma
!!!  9-function_eval_sigma_v1             !! computes sgma with some adaptivity (constant alpha)
!!!                                       !! and using brute force
!!!
!!! INDEX OF DERIVED TYPES:
!!!
!!! Public:
!!!
!!! 1-Fast_Sigma_stuff       !! Contains all the information to evaluate sgma, including the tree and some constants
!!!
!!! Private: (none)
!!!
!!! MODULES REQUIRED:
!!!
!!! 1-Mod_TreeLRD
!!! 2-ModType_Smooth_Surface
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module Mod_Fast_Sigma
!! This subroutine evals the fast gaussian on a vector of targets
use Mod_TreeLRD
use ModType_Smooth_Surface
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!    DEFINITION OF TYPES     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type, public :: Fast_Sigma_stuff
    type ( TreeLRD ), pointer :: TreeLRD_1 => null ()
    real ( kind = 8 ) :: sgma_mean                          !! mean value of the vector W_sgmas_mem
    real ( kind = 8 ) :: alpha                          !! mean value of the vector W_sgmas_mem
end type Fast_Sigma_stuff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!    END DEFINITION OF TYPES     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!    DEFINITION OF SUBROUTINES     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


public :: setup_tree_sigma_geometry
public :: function_eval_sigma

private :: find_centroids_sigmas
private :: fast_gaussian_global
private :: fast_gaussian_box
private :: fast_gaussian_box_grad
public :: initial_guess_sgma
public :: eval_sigma_erf
private :: function_eval_sigma_v1
private :: function_eval_sigma_v0

contains





  
subroutine setup_tree_sigma_geometry(FSS_1,Geometry1,rlam)
implicit none
!!This routine initialize the tree to allow a fast sigma evaluator

    !List of calling arguments
    type ( Geometry ), intent(inout) :: Geometry1       !! Geometry information
    type ( Fast_Sigma_stuff ), pointer ::  FSS_1             !! Tree where we place the centroids and sigma values at the centroids

    !List of local variables
    integer n_max_leaf
    double precision :: Centroids(3,Geometry1%ntri)
    double precision :: sgmas(Geometry1%ntri), rlam, sigma0

    !! Tree where we place the centroids and sigma values at the
    !! centroids
    type ( TreeLRD ), pointer ::  TreeLRD_1             
    integer n_leaf_boxes,count1
    real ( kind = 8 ), allocatable :: centers(:,:),sgmas_2(:)

        allocate(FSS_1)
        n_max_leaf=1
        call find_centroids_sigmas(Geometry1,Centroids,sgmas,rlam)   !! Find the centroids and sgmas on each centroid
!!!! ESTO ESTÁ CAMBIADO adapt_flag=1
        call start_tree(FSS_1%TreeLRD_1,n_max_leaf,Centroids,sgmas,Geometry1%ntri_sk) !! Allocate centroids and sgmas on the tree
!!!! ESTO ESTÁ CAMBIADO adapt_flag=2
!call start_tree(FSS_1%TreeLRD_1,n_max_leaf,Geometry1%skeleton_Points,Geometry1%skeleton_Points(1,:),Geometry1%n_Sk_points) !! Allocate centroids and sgmas on the tree

        call make_level_restricted(FSS_1%TreeLRD_1)                   !! Make tree level restricted
        call setup_nearest(FSS_1%TreeLRD_1%Main_box)
        call defrag_tree_Points(FSS_1%TreeLRD_1)                      !! Defragment the tree
        !write (*,*) 'Start P_vect setup'
!        write (*,*) 'big box data (center): ', FSS_1%TreeLRD_1%Main_Box%Box_center
!        write (*,*) 'big box data (size): ', FSS_1%TreeLRD_1%Main_Box%Box_size

!write (*,*) 'x_min: ',FSS_1%TreeLRD_1%Main_Box%Box_center(1)-FSS_1%TreeLRD_1%Main_Box%Box_size/2.0d0
!write (*,*) 'x_max: ',FSS_1%TreeLRD_1%Main_Box%Box_center(1)+FSS_1%TreeLRD_1%Main_Box%Box_size/2.0d0
!write (*,*) 'y_min: ',FSS_1%TreeLRD_1%Main_Box%Box_center(2)-FSS_1%TreeLRD_1%Main_Box%Box_size/2.0d0
!write (*,*) 'y_max: ',FSS_1%TreeLRD_1%Main_Box%Box_center(2)+FSS_1%TreeLRD_1%Main_Box%Box_size/2.0d0
!write (*,*) 'z_min: ',FSS_1%TreeLRD_1%Main_Box%Box_center(3)-FSS_1%TreeLRD_1%Main_Box%Box_size/2.0d0
!write (*,*) 'z_max: ',FSS_1%TreeLRD_1%Main_Box%Box_center(3)+FSS_1%TreeLRD_1%Main_Box%Box_size/2.0d0


        call setup_P_vect(FSS_1%TreeLRD_1,FSS_1%TreeLRD_1%Main_Box%Box_center,FSS_1%TreeLRD_1%Main_Box%Box_size)

        !write (*,*) 'finito P_vect setup'
!        read (*,*)

        !Compute the mean value of the sigmas to allow nonadaptivity
        FSS_1%sgma_mean = sum(FSS_1%TreeLRD_1%W_sgmas_mem)/Geometry1%ntri
        !FSS_1%alpha = 1.0d0/(2.0d0*maxval(FSS_1%TreeLRD_1%W_sgmas_mem)**2)
        !FSS_1%alpha=1.0d0/maxval(FSS_1%TreeLRD_1%W_sgmas_mem)**2/5.0d0
        !print *, 'FSS alpha = ', FSS_1%alpha
        
        sigma0 = sqrt(5/2.0d0)*maxval(FSS_1%TreeLRD_1%W_sgmas_mem)
        !print *, 'sigma0 = ', sigma0

        FSS_1%alpha = 1/sigma0**2/2
        !print *, 'using sigma0, FSS alpha = ', FSS_1%alpha
        !stop

!        n_leaf_boxes=total_number_leaf_boxes(FSS_1%TreeLRD_1%Main_box)
!        allocate(centers(3,n_leaf_boxes))
!        allocate(sgmas_2(n_leaf_boxes))
!        call all_centers_sgmas(FSS_1%TreeLRD_1,centers,sgmas_2,n_leaf_boxes)
!        call deallocate_tree(FSS_1%TreeLRD_1)
!        call start_tree(FSS_1%TreeLRD_1,n_max_leaf,centers,sgmas_2,n_leaf_boxes) !! Allocate centroids and sgmas on the tree
!        call make_level_restricted(FSS_1%TreeLRD_1)                   !! Make tree level restricted
!        call setup_nearest(FSS_1%TreeLRD_1%Main_box)
!        call defrag_tree_Points(FSS_1%TreeLRD_1)                      !! Defragment the tree
!        deallocate(centers)
!        deallocate(sgmas_2)
!        FSS_1%sgma_mean=sum(FSS_1%TreeLRD_1%W_sgmas_mem)/Geometry1%ntri
!        FSS_1%alpha=1.0d0/maxval(FSS_1%TreeLRD_1%W_sgmas_mem)**2/5.0d0
!        FSS_1%alpha=FSS_1%alpha/4.0d0

 !       write (*,*) 'alpha val: ',FSS_1%alpha,1/dsqrt(2*FSS_1%alpha),20.0d0/sqrt(2.0d0*FSS_1%alpha)


return
end






subroutine function_eval_sigma(FSS_1,targ_vect,n_targ,sgma, &
    sgma_x,sgma_y,sgma_z,adapt_flag)
  implicit none
  
  !! Use adapt_flag=0 for completely non-adaptive case (sgma=constant
  !! across the geometry) (all derivatives zero)

  !! Use adapt_flag=1 for some adaptivity (constant alpha) it will be
  !! faster in case of no need of HUGE adaptivity

  !! Use adapt_flag=2 for HUGE adaptivity requirement

  !List of calling arguments

  !! This is the Tree generated by setup_tree_sigma_geometry containing
  type ( Fast_Sigma_stuff ), pointer ::  FSS_1

  !! points of centroids and sigmas
  !! Number of targets and flags to specify adaptivity and speed
  !! adapt_flag=0 means constant sigma (no adaptivity at all)
  !! adapt_flag=1 means constant alfa (adaptiviti)
  !!
  !! adapt_flag=2 means recursive sigma definition (HUGE adaptivity
  !!   and a bit slower)
  integer, intent(in) :: n_targ,adapt_flag

  !! vector with the targets
  real ( kind = 8 ), intent(in) :: targ_vect(3,n_targ)

  !! output sigma and gradient of sigma
  real ( kind = 8 ), intent(out) :: sgma(n_targ),sgma_x(n_targ)
  double precision :: sgma_y(n_targ),sgma_z(n_targ) 

  !! this function calls the corresponding subroutine with the
  !! appropriate sigma generator for the selected flats

  !
  ! sigma will be constant across the geometry
  !
  if (adapt_flag .eq. 0) then
    call function_eval_sigma_v0(targ_vect,n_targ,FSS_1,sgma, &
        sgma_x,sgma_y,sgma_z)
    return
  endif

  !
  ! otherwise sigma is adaptive
  !
  if ( (adapt_flag .eq. 1) .or. (adapt_flag .eq. 2) ) then
    call fast_gaussian_global_new(FSS_1,targ_vect,n_targ,sgma, &
        sgma_x,sgma_y,sgma_z, adapt_flag)
    return
  endif

  !
  ! otherwise flag error
  !
  print *, 'stopping, adapt_flag = ', adapt_flag
  stop

  return
end subroutine function_eval_sigma






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine find_centroids_sigmas(Geometry1,Centroids,sgmas,rlam)
implicit none

    !List of calling arguments
    type (Geometry), intent(inout) :: Geometry1
    real ( kind = 8 ), intent(inout) ::  Centroids(3,Geometry1%ntri),sgmas(Geometry1%ntri)


    !List of local variables
    integer count
    real ( kind = 8 ) lambda,d1,d2,d3,dmax,sgma_mean,rlam
    real ( kind = 8 ) P1(3),P2(3),P3(3),Pc(3),Pd1(3),Pd2(3),Pd3(3)


        sgma_mean=0.0d0
        do count=1,Geometry1%ntri
            P1=Geometry1%Points(:,Geometry1%Tri(1,count))
            P2=Geometry1%Points(:,Geometry1%Tri(2,count))
            P3=Geometry1%Points(:,Geometry1%Tri(3,count))
            Pc=(P1+P2+P3)/3.0d0
            Pd1=P2-P1
            Pd2=P3-P2
            Pd3=P3-P1
            d1=sqrt(Pd1(1)**2+Pd1(2)**2+Pd1(3)**2)
            d2=sqrt(Pd2(1)**2+Pd2(2)**2+Pd2(3)**2)
            d3=sqrt(Pd3(1)**2+Pd3(2)**2+Pd3(3)**2)
            dmax=max(d1,d2,d3)
            Centroids(:,count)=Pc
            sgmas(count)=dmax/rlam
        enddo
!!        if (adapt_flag==0) then
!!            write (*,*) 'adaptive flag = 0'
!!            sgma_mean=sum(Geometry1%sgmas)/Geometry1%ntri
!!            do count=1,Geometry1%ntri
!!                Geometry1%sgmas(count)=sgma_mean
!!            enddo
!!        endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine function_eval_sigma_v0(targ_vect,n_targ,FSS_1,sgma,Gx,Gy,Gz)
implicit none
!! This subroutine computes sgma with no adaptivity at all (constant sgma and derivatives=0)

    !List of calling arguments
    type ( Fast_Sigma_stuff ), pointer ::  FSS_1
    integer, intent(in) :: n_targ
    real ( kind = 8 ), intent(in) :: targ_vect(3,n_targ)
    real ( kind = 8 ), intent(out) :: sgma(n_targ),Gx(n_targ),Gy(n_targ),Gz(n_targ)

    !List of local variables
    real ( kind = 8 ) F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z,di,my_exp
    integer count1,count2

        do count1=1,n_targ
            Gx(count1)=0.0d0
            Gy(count1)=0.0d0
            Gz(count1)=0.0d0
            sgma(count1)=FSS_1%sgma_mean
        enddo

return
end subroutine function_eval_sigma_v0







subroutine fast_gaussian_global_new(FSS_1, targ_vect, n_targ, sgma, &
    sgma_x, sgma_y, sgma_z, adapt_flag)

  implicit none

  !
  !! This program runs the fixed point iteration to find sigma
  !! starting from initial guess that comes from a the size of the box
  !! in which the target would land in the level restricted tree.
  !! Notice that the tree covers the full domain
  !
  
  ! List of calling arguments
  type ( Fast_Sigma_stuff ), pointer ::  FSS_1

  integer :: n_targ, adapt_flag
  double precision :: targ_vect(3,n_targ)
  double precision :: sgma(n_targ), sgma_x(n_targ)
  double precision :: sgma_y(n_targ),sgma_z(n_targ)

  ! List of local variables
  integer :: count1, count2
  double precision :: sgma_min, sgma_max
  double precision :: d_aux,alpha,my_exp,F,D,dF_x,dF_y,dF_z
  double precision :: dF_sigma,dD_sigma,denom
  double precision :: dD_x,dD_y,dD_z,sgm_rad,tol,pot1,pot2,err
  double precision :: alpha0, alpham1
  double precision :: sigma0, sigmam1
  double precision :: sgm_rad0,sgm_radm1
  double precision :: f0,fm1,d0,dm1,dd

  integer nhalf,j
  
  double precision :: fprime, tes,tes2,resid,hstep,sigma

  integer :: count_max
  type ( TreeLRD ), pointer :: TreeLRD_1

  TreeLRD_1 => FSS_1%TreeLRD_1

  sgma_min = minval(TreeLRD_1%W_sgmas_mem)
  sgma_max = maxval(TreeLRD_1%W_sgmas_mem)

!  call prin2('sgma_min=*',sgma_min,1)
!  call prin2('sgma_max=*',sgma_max,1)

  tol = 1.0d-14
  nhalf = 3

  count_max = (log(sgma_max-sgma_min)-log(tol))/log(2.0d0)+4

!  call prinf('count_max=*',count_max,1)
  


  !!!!$ call omp_set_num_threads(1)

  !$OMP PARALLEL DO DEFAULT(SHARED), &
  !$OMP PRIVATE(count1,count2,err,sgm_rad,alpha,F,D), &
  !$OMP PRIVATE(dF_x,dF_y,dF_z,dD_x,dD_y,dD_z,pot1,pot2), &
  !$OMP PRIVATE(alpha0,alpham1,sgm_rad0,sgm_radm1,sigma), &
  !$OMP PRIVATE(sigmam1,sigma0,resid,tes,tes2,fprime), &
  !$OMP PRIVATE(f0,d0,fm1,dm1,hstep,j,dF_sigma,dD_sigma,denom,dd)
  do count1 = 1,n_targ

    !! Esto hay que revisarlo!!!!
    !! This has to be reviewed!! !
    alpha = FSS_1%alpha
    sgm_rad = 12.0d0/sqrt(2.0d0*alpha)

    !call prin2('. . . before adapt 2*', alpha, 0)
    !call prin2('sgm_rad = *', sgm_rad, 1)
    !call prin2('alpha = *', alpha, 1)

    
    if (adapt_flag .eq. 2) then


 !
 !!     Get two initial guesses for secant code
 !

      err=tol+1.0d0

      alpha0 = alpha
      alpham1 = 2.0d0*alpha

      sgm_rad0 = 12.0d0/sqrt(2.0d0*alpha0)
      sgm_radm1 = 12.0d0/sqrt(2.0d0*alpham1)

      call fast_gaussian_box_v2(TreeLRD_1%Main_box, &
            targ_vect(:,count1), sgm_rad0, alpha0, f0, d0)
      call fast_gaussian_box_v2(TreeLRD_1%Main_box, &
            targ_vect(:,count1), sgm_radm1, alpham1, fm1, dm1)


!
!!    sigma is current guess, sigmam1 is previous guess
!

      sigma = f0/d0
      sigmam1 = fm1/dm1

      alpham1 = 1/(2.0d0*sigmam1**2)
      sgm_radm1 = 12.0d0/sqrt(2.0d0*alpham1)
      call fast_gaussian_box_v2(TreeLRD_1%Main_box, &
            targ_vect(:,count1), sgm_radm1, alpham1, fm1, dm1)




      !!!! call initial_guess_sgma(TreeLRD_1%Main_box,targ_vect(:,count1),sgm_rad)
      !!!! alpha=1/(10*sgm_rad/70.71d0)**2


      !call prin2('. . . inside adapt 2*', alpha, 0)
      count2=0

      do count2 = 1,count_max
        sigma0 = sigma

        alpha0 = 1/(2.0d0*sigma0**2)
        sgm_rad0 = 12.0d0/sqrt(2.0d0*alpha0)
        call fast_gaussian_box_v2(TreeLRD_1%Main_box, &
            targ_vect(:,count1), sgm_rad0, alpha0, f0, d0)
        
        resid = f0/d0 - sigma0

        tes = abs(resid)


        if(tes.le.tol) goto 1001

        fprime = (f0/d0-fm1/dm1)/(sigma0-sigmam1)

        fm1 = f0
        dm1 = d0
        
        sigmam1 = sigma0

        hstep = -resid/(-1.0d0+fprime)
        do j=1,nhalf
           sigma = sigma0 + hstep
           alpha = 1.0d0/(2.0d0*sigma**2)
           sgm_rad = 12.0d0/sqrt(2.0d0*alpha)
           call fast_gaussian_box_v2(TreeLRD_1%Main_box, &
              targ_vect(:,count1), sgm_rad, alpha, f, d)
           
           resid = f/d-sigma
           tes2 = abs(resid)
           if(tes2.le.tol) goto 1001
           if(tes2.lt.tes) then
             goto 1000
           else
             hstep = hstep*0.5d0
           endif
        enddo
        print *, "Exceeding max damping steps in sigma evaluator" 
 1000 continue

      end do
 1001 continue      

      
      alpha = 1.0d0/(2*sigma**2)
      sgm_rad = 12.0d0/sqrt(2.0d0*alpha)
      
    endif

    !
    ! using the obtained value, call the fast Gauss transform
    !
  
    call fast_gaussian_box_grad(adapt_flag,TreeLRD_1%Main_box, &
        targ_vect(:,count1),sgm_rad,alpha,F,D, &
        dF_x,dF_y,dF_z,dD_x,dD_y,dD_z, dF_sigma, dD_sigma)

    sgma(count1) = F/D
    if(adapt_flag.eq.2) then
      if(abs(sigma-sgma(count1)).ge.tol) then
        print *, "sigma did not converge for point=",count1
        stop
      endif
    endif

    sgma_x(count1) = (dF_x*D-F*dD_x)/D**2
    sgma_y(count1) = (dF_y*D-F*dD_y)/D**2
    sgma_z(count1) = (dF_z*D-F*dD_z)/D**2

    if(adapt_flag.eq.2) then

      dd = (D*dF_sigma - F*dD_sigma)/D**2/sgma(count1)**3
      denom = 1.0d0/(1.0d0-dd)

      sgma_x(count1) = sgma_x(count1)*denom
      sgma_y(count1) = sgma_y(count1)*denom
      sgma_z(count1) = sgma_z(count1)*denom
    endif
  enddo
  !$omp end parallel do


  return
end subroutine fast_gaussian_global_new





subroutine fast_gaussian_global(FSS_1, targ_vect, n_targ, sgma, &
    sgma_x, sgma_y, sgma_z, adapt_flag)

  implicit none

  !
  !! This program runs the fixed point iteration to find sigma
  !! starting from initial guess that comes from a the size of the box
  !! in which the target would land in the level restricted tree.
  !! Notice that the tree covers the full domain
  !
  
  ! List of calling arguments
  type ( Fast_Sigma_stuff ), pointer ::  FSS_1

  integer, intent(in) :: n_targ,adapt_flag
  real ( kind = 8 ), intent(in) :: targ_vect(3,n_targ)
  real ( kind = 8 ), intent(out) :: sgma(n_targ),sgma_x(n_targ)
  double precision :: sgma_y(n_targ),sgma_z(n_targ)

  ! List of local variables
  integer count1,count2
  double precision :: d_aux,alpha,my_exp,F,D,dF_x,dF_y,dF_z
  double precision :: dD_x,dD_y,dD_z,sgm_rad,tol,pot1,pot2,err
  double precision :: dF_sigma, dD_sigma,denom

  type ( TreeLRD ), pointer :: TreeLRD_1

  TreeLRD_1 => FSS_1%TreeLRD_1

  tol = 1.0d-10

  !$OMP PARALLEL DO DEFAULT(SHARED), &
  !$OMP& PRIVATE(count1,count2,err,sgm_rad,alpha,F,D), &
  !$OMP& PRIVATE(dF_x,dF_y,dF_z,dD_x,dD_y,dD_z,pot1,pot2,dF_sigma), &
  !$OMP& PRIVATE(dD_sigma,denom)
  do count1 = 1,n_targ

    if (adapt_flag==2) then
      err=tol+1.0d0
      !                call initial_guess_sgma(TreeLRD_1%Main_box,targ_vect(:,count1),sgm_rad)
      !alpha=1/(10*sgm_rad/70.71d0)**2

      !! Esto hay que revisarlo!!!!
      alpha = FSS_1%alpha
      sgm_rad = 12.0d0/sqrt(2.0d0*alpha)

      count2=0
      do while ( (err>tol) .and. (count2<1) )
        call fast_gaussian_box_v2(TreeLRD_1%Main_box, &
            targ_vect(:,count1), sgm_rad, alpha, F, D)

        pot1=F/D
        ! write (*,*) 'pot iteration: ', pot1,pot2,err

        err=dabs(pot2-pot1)/dabs(pot1)
        count2 = count2+1
        pot2=pot1
        sgm_rad = 70.7d0*pot2
        alpha = 1/(10*pot2)**2
      end do

    else
      !! Esto hay que revisarlo!! !!
      !! This has to be reviewed!! !
      
      !print *, 'adapt_flag = ', adapt_flag
      !stop
      alpha = FSS_1%alpha
      sgm_rad = 12.0d0/sqrt(2.0d0*alpha)
    endif

    call fast_gaussian_box_grad(adapt_flag, TreeLRD_1%Main_box,targ_vect(:,count1),&
      sgm_rad,alpha,F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z,dF_sigma, dD_sigma)

    sgma(count1) = F/D
    sgma_x(count1) = (dF_x*D-F*dD_x)/D**2
    sgma_y(count1) = (dF_y*D-F*dD_y)/D**2
    sgma_z(count1) = (dF_z*D-F*dD_z)/D**2

    if(adapt_flag.eq.2) then
      denom = 1.0d0/(1.0d0-(D*dF_sigma - F*dD_sigma)/D**2/sgma(count1)**3)
      sgma_x(count1) = sgma_x(count1)*denom
      sgma_y(count1) = sgma_y(count1)*denom
      sgma_z(count1) = sgma_z(count1)*denom
    endif
  enddo
  !$OMP END PARALLEL DO


  return
end subroutine fast_gaussian_global







recursive subroutine fast_gaussian_box(Current_box,targ,sgm_rad,alpha,F,D)
implicit none
!! This subroutine evals the fast gaussian recursively for one target (F and D)

!! Main_box: (inout) pointer to null where the main box will be allocated
!! n_max_leaf: (in) maximum number of points allowed in a box to stop splitting
!! targ: (in) (3) real vector that contains the target
!! sgm_rad: (in) radius associated with the target
!! F, D: (inout) numerator and denominator of the solution

    !List of calling arguments
    type ( Box ), pointer :: Current_box
    real ( kind = 8 ), intent(in) :: targ(3),sgm_rad,alpha
    real ( kind = 8 ), intent(inout) :: F,D

    !List of local variables
    integer :: count1,count2,count3
    real ( kind = 8 ) d_aux,my_exp
    real ( kind = 8 ) d2

        if (.not. (associated(Current_box%Parent))) then
            F=0.0d0
            D=0.0d0
        endif

        if (((Current_box%Box_center(1)+Current_box%Box_size/2.0d0) >= (targ(1)-sgm_rad)).and. &
            & ((Current_box%Box_center(1)-Current_box%Box_size/2.0d0) <= (targ(1)+sgm_rad)) .and. &
            ((Current_box%Box_center(2)+Current_box%Box_size/2.0d0) >= (targ(2)-sgm_rad)).and. &
            & ((Current_box%Box_center(2)-Current_box%Box_size/2.0d0) < (targ(2)+sgm_rad)) .and. &
            ((Current_box%Box_center(3)+Current_box%Box_size/2.0d0) >= (targ(3)-sgm_rad)).and. &
            & ((Current_box%Box_center(3)-Current_box%Box_size/2.0d0) < (targ(3)+sgm_rad))) then

            if ((Current_box%n_points)==0) then
                do count1=1,8
                    if (associated(Current_box%Children(count1)%BP)) then
                        call fast_gaussian_box(Current_box%Children(count1)%BP,targ,sgm_rad,alpha,F,D)
                    endif
                enddo
            else
!               alpha=1.0d0/(2*(sgm_rad/10.0d0))
                do count1=1,Current_box%n_points
                    d2=(Current_box%Points(1,count1)-targ(1))**2+(Current_box%Points(2,count1)&
                     &-targ(2))**2+(Current_box%Points(3,count1)-targ(3))**2
                    my_exp=exp(-alpha*d2)
                    F=F+Current_box%sgmas(count1)*my_exp
                    D=D+my_exp
                enddo
            endif
        endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine fast_gaussian_box_v2(Current_box,targ,sgm_rad,alpha,F,D)
  implicit none
  
  !! This subroutine evals the fast gaussian recursively for one
  !! target (F and D)

  !! Current_box: (inout) pointer to null where the main box will be allocated
  !! n_max_leaf: (in) maximum number of points allowed in a box to stop splitting
  !! targ: (in) (3) real vector that contains the target
  !! sgm_rad: (in) radius associated with the target
  !! F, D: (inout) numerator and denominator of the solution
  !! dF_x,dF_y,dF_z,dD_x,dD_y,dD_z all the derivatives

  !List of calling arguments
  type ( Box ), pointer :: Current_box
  real ( kind = 8 ), intent(in) :: targ(3),sgm_rad,alpha
  real ( kind = 8 ), intent(inout) :: F,D

    !List of local variables
    integer  count1,count2,count3
    real ( kind = 8 ) d_aux,my_exp
    real ( kind = 8 ) d2

        if (.not. (associated(Current_box%Parent))) then
            F=0.0d0
            D=0.0d0
        endif

        if (((Current_box%Box_center(1)+Current_box%Box_size/2.0d0) >= (targ(1)-sgm_rad)).and. &
            & ((Current_box%Box_center(1)-Current_box%Box_size/2.0d0) <= (targ(1)+sgm_rad)) .and. &
            ((Current_box%Box_center(2)+Current_box%Box_size/2.0d0) >= (targ(2)-sgm_rad)).and. &
            & ((Current_box%Box_center(2)-Current_box%Box_size/2.0d0) < (targ(2)+sgm_rad)) .and. &
            ((Current_box%Box_center(3)+Current_box%Box_size/2.0d0) >= (targ(3)-sgm_rad)).and. &
            & ((Current_box%Box_center(3)-Current_box%Box_size/2.0d0) < (targ(3)+sgm_rad))) then
            if ((Current_box%n_points)==0) then
                do count1=1,8
                    if (associated(Current_box%Children(count1)%BP)) then
                        call fast_gaussian_box_v2(Current_box%Children(count1)%BP,targ,sgm_rad,alpha,F,D)
                    endif
                enddo
            else
                do count1=1,Current_box%n_points
                    d2=(Current_box%Points(1,count1)-targ(1))**2+(Current_box%Points(2,count1)&
                     &-targ(2))**2+(Current_box%Points(3,count1)-targ(3))**2
                    my_exp=exp(-alpha*d2)
                    F=F+Current_box%sgmas(count1)*my_exp
                    D=D+my_exp
                enddo
            endif
        endif

return
end




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine fast_gaussian_box_grad(adapt_flag,Current_box,targ,sgm_rad,&
    alpha,F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z,dF_sigma, dD_sigma)
  implicit none
  !
  !! This subroutine evals the fast gaussian recursively for one
  !! target (F and D)

  !! Current_box: (inout) pointer to null where the main box will be
  !! allocated n_max_leaf: (in) maximum number of points allowed in a
  !! box to stop splitting targ: (in) (3) real vector that contains
  !! the target sgm_rad: (in) radius associated with the target F, D:
  !! (inout) numerator and denominator of the solution
  !! dF_x,dF_y,dF_z,dD_x,dD_y,dD_z all the derivatives
  !
  !! dF_sigma,dD_sigma, are the sigma derivatives.
  !! Note, we are missing a factor of sigma**3 in dF_sigma and dD_sigma
  !!

  !List of calling arguments
  type ( Box ), pointer :: Current_box
  real ( kind = 8 ), intent(in) :: targ(3),sgm_rad,alpha
  real ( kind = 8 ), intent(inout) :: F,D,dF_x,dF_y,dF_z
  double precision :: dD_x,dD_y,dD_z,dF_sigma,dD_sigma

  !List of local variables
  integer  count1,count2,count3, adapt_flag
  real ( kind = 8 ) d_aux,my_exp
  real ( kind = 8 ) d2, rtmp

  if (.not. (associated(Current_box%Parent))) then
    F=0.0d0
    D=0.0d0
    dF_x=0.0d0
    dF_y=0.0d0
    dF_z=0.0d0
    dD_x=0.0d0
    dD_y=0.0d0
    dD_z=0.0d0
    dF_sigma = 0.0d0
    dD_sigma = 0.0d0
  endif

  if (((Current_box%Box_center(1)+Current_box%Box_size/2.0d0) >= (targ(1)-sgm_rad)).and. &
       ((Current_box%Box_center(1)-Current_box%Box_size/2.0d0) <= (targ(1)+sgm_rad)) .and. &
      ((Current_box%Box_center(2)+Current_box%Box_size/2.0d0) >= (targ(2)-sgm_rad)).and. &
       ((Current_box%Box_center(2)-Current_box%Box_size/2.0d0) < (targ(2)+sgm_rad)) .and. &
      ((Current_box%Box_center(3)+Current_box%Box_size/2.0d0) >= (targ(3)-sgm_rad)).and. &
       ((Current_box%Box_center(3)-Current_box%Box_size/2.0d0) < (targ(3)+sgm_rad))) then

    ! inside of a box
    
    if ((Current_box%n_points)==0) then
      do count1=1,8
        if (associated(Current_box%Children(count1)%BP)) then
          call fast_gaussian_box_grad(adapt_flag,&
              Current_box%Children(count1)%BP,targ,sgm_rad,alpha,&
              F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z,dF_sigma,dD_sigma)
        endif
      enddo
    else

      ! alpha=1.0d0/(2*(sgm_rad/10.0d0))

      do count1=1,Current_box%n_points
        d2=(Current_box%Points(1,count1)-targ(1))**2 &
            +(Current_box%Points(2,count1) &
            -targ(2))**2+(Current_box%Points(3,count1)-targ(3))**2
        my_exp = exp(-alpha*d2)
        rtmp = my_exp*Current_box%sgmas(count1)


        F=F+rtmp
        D=D+my_exp
        dF_x=dF_x+(-2*alpha)*rtmp*(targ(1)-Current_box%Points(1,count1))
        dF_y=dF_y+(-2*alpha)*rtmp*(targ(2)-Current_box%Points(2,count1))
        dF_z=dF_z+(-2*alpha)*rtmp*(targ(3)-Current_box%Points(3,count1))
        dD_x=dD_x+(-2*alpha)*my_exp*(targ(1)-Current_box%Points(1,count1))
        dD_y=dD_y+(-2*alpha)*my_exp*(targ(2)-Current_box%Points(2,count1))
        dD_z=dD_z+(-2*alpha)*my_exp*(targ(3)-Current_box%Points(3,count1))

        if(adapt_flag.eq.2) then
          dF_sigma = dF_sigma + d2*rtmp
          dD_sigma = dD_sigma + d2*my_exp
        endif
      enddo
    endif
  else

    ! out of a box

  endif

  return
end subroutine fast_gaussian_box_grad





recursive subroutine initial_guess_sgma(Current_box,targ,sgm_rad)
implicit none
!! This function provides the sizee of the leaf box where a target point (non-existing on the tree)
!! would land

    !List of calling arguments
    type ( Box ), pointer :: Current_box
    real ( kind = 8 ), intent(in) :: targ(3)
    real ( kind = 8 ), intent(out) :: sgm_rad

    !List of local variables

        if (Current_box%is_leaf) then
            sgm_rad=Current_box%Box_size*1.5d0/5.0d0
        else
            if (targ(1)>Current_box%Box_center(1)) then
                if (targ(2)>Current_box%Box_center(2)) then
                    if (targ(3)>Current_box%Box_center(3)) then
                        call initial_guess_sgma(Current_box%Children(1)%BP,targ,sgm_rad)
                    else
                        call initial_guess_sgma(Current_box%Children(5)%BP,targ,sgm_rad)
                    endif
                else
                    if (targ(3)>Current_box%Box_center(3)) then
                        call initial_guess_sgma(Current_box%Children(3)%BP,targ,sgm_rad)
                    else
                        call initial_guess_sgma(Current_box%Children(7)%BP,targ,sgm_rad)
                    endif
                endif
            else
                if (targ(2)>Current_box%Box_center(2)) then
                    if (targ(3)>Current_box%Box_center(3)) then
                        call initial_guess_sgma(Current_box%Children(2)%BP,targ,sgm_rad)
                    else
                        call initial_guess_sgma(Current_box%Children(6)%BP,targ,sgm_rad)
                    endif
                else
                    if (targ(3)>Current_box%Box_center(3)) then
                        call initial_guess_sgma(Current_box%Children(4)%BP,targ,sgm_rad)
                    else
                        call initial_guess_sgma(Current_box%Children(8)%BP,targ,sgm_rad)
                    endif
                endif
            endif
        endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine function_eval_sigma_v1(targ_vect,n_targ,FSS_1,sgma,Gx,Gy,Gz)
implicit none
!! This subroutine computes sgma by brute force with some adaptivity (constant alpha)


    !List of calling arguments
    type ( Fast_Sigma_stuff ), pointer ::  FSS_1
    integer, intent(in) :: n_targ
    real ( kind = 8 ), intent(in) :: targ_vect(3,n_targ)
    real ( kind = 8 ), intent(out) :: sgma(n_targ),Gx(n_targ),Gy(n_targ),Gz(n_targ)

    !List of local variables
    real ( kind = 8 ) F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z,di,my_exp,alpha
    integer count1,count2
    type ( TreeLRD ), pointer :: TreeLRD_1


        TreeLRD_1 => FSS_1%TreeLRD_1
!        alpha=1.0d0/maxval(TreeLRD_1%W_sgmas_mem)**2/5.0d0
        alpha=FSS_1%alpha
        write (*,*) 'doing it slow'
        do count1=1,n_targ
            F=0.0d0
            D=0.0d0
            dF_x=0.0d0
            dF_y=0.0d0
            dF_z=0.0d0
            dD_x=0.0d0
            dD_y=0.0d0
            dD_z=0.0d0
            do count2=1,TreeLRD_1%ntot_W_Pts_mem
                di=sqrt((TreeLRD_1%W_Pts_mem(1,count2)-targ_vect(1,count1))**2+(TreeLRD_1%W_Pts_mem(2,count2)&
                 &-targ_vect(2,count1))**2+(TreeLRD_1%W_Pts_mem(3,count2)-targ_vect(3,count1))**2)
                my_exp=exp(-alpha*di**2)
                F=F+TreeLRD_1%W_sgmas_mem(count2)*my_exp
                dF_x=dF_x+(-2*alpha)*TreeLRD_1%W_sgmas_mem(count2)*my_exp*(targ_vect(1,count1)-TreeLRD_1%W_Pts_mem(1,count2))
                dF_y=dF_y+(-2*alpha)*TreeLRD_1%W_sgmas_mem(count2)*my_exp*(targ_vect(2,count1)-TreeLRD_1%W_Pts_mem(2,count2))
                dF_z=dF_z+(-2*alpha)*TreeLRD_1%W_sgmas_mem(count2)*my_exp*(targ_vect(3,count1)-TreeLRD_1%W_Pts_mem(3,count2))
                dD_x=dD_x+(-2*alpha)*my_exp*(targ_vect(1,count1)-TreeLRD_1%W_Pts_mem(1,count2))
                dD_y=dD_y+(-2*alpha)*my_exp*(targ_vect(2,count1)-TreeLRD_1%W_Pts_mem(2,count2))
                dD_z=dD_z+(-2*alpha)*my_exp*(targ_vect(3,count1)-TreeLRD_1%W_Pts_mem(3,count2))
                D=D+my_exp
            enddo
            Gx(count1)=(dF_x*D-F*dD_x)/D**2
            Gy(count1)=(dF_y*D-F*dD_y)/D**2
            Gz(count1)=(dF_z*D-F*dD_z)/D**2
            sgma(count1)=F/D
!            write (*,*) 'F/D slow: ',F,D
!            read (*,*)
        enddo

return
end



subroutine eval_sigma_erf(Points,n,FSS_1,sgma,sgma_x,sgma_y,sgma_z)
implicit none

    !List of calling arguments
    type ( Fast_Sigma_stuff ), pointer ::  FSS_1
!    type ( TreeLRD ), pointer :: TreeLRD_1
    integer, intent(in) :: n
    real ( kind = 8 ), intent(in) :: Points(3,n)
    real ( kind = 8 ), intent(out) :: sgma(n),sgma_x(n),sgma_y(n),sgma_z(n)

    !List of local variables
    integer count1
    type ( Box ), pointer :: p_box
    real ( kind = 8 ) x,y,z

    do count1=1,n
        call find_leaf_node(FSS_1%TreeLRD_1%Main_box,Points(:,count1),p_box)
        if (p_box%P_vect(2,2,2)<0.0d0) then
            x=(Points(1,count1)-p_box%Box_center(1))/p_box%Box_size+0.5d0
            y=(Points(2,count1)-p_box%Box_center(2))/p_box%Box_size+0.5d0
            z=(Points(3,count1)-p_box%Box_center(3))/p_box%Box_size+0.5d0
            sgma(count1)=interp_3D(x,y,z,p_box%P_vect)
            call interp_3D_grad(x,y,z,p_box%P_vect,sgma(count1),sgma_x(count1),sgma_y(count1),sgma_z(count1))
        else
            sgma(count1)=p_box%P_vect(2,2,2)
            sgma_x(count1)=0.0d0
            sgma_y(count1)=0.0d0
            sgma_z(count1)=0.0d0
        endif
    enddo

return
end


recursive subroutine find_leaf_node(Current_box,targ,p_box)
implicit none
!! This function provides the size of the leaf box where a target point (non-existing on the tree)
!! would land

    !List of calling arguments
    type ( Box ), pointer :: Current_box
    type ( Box ), pointer :: p_box
    real ( kind = 8 ), intent(in) :: targ(3)

    !List of local variables

        if (Current_box%is_leaf) then
            p_box => Current_box
        else
            if (targ(1)>Current_box%Box_center(1)) then
                if (targ(2)>Current_box%Box_center(2)) then
                    if (targ(3)>Current_box%Box_center(3)) then
                        call find_leaf_node(Current_box%Children(1)%BP,targ,p_box)
                    else
                        call find_leaf_node(Current_box%Children(5)%BP,targ,p_box)
                    endif
                else
                    if (targ(3)>Current_box%Box_center(3)) then
                        call find_leaf_node(Current_box%Children(3)%BP,targ,p_box)
                    else
                        call find_leaf_node(Current_box%Children(7)%BP,targ,p_box)
                    endif
                endif
            else
                if (targ(2)>Current_box%Box_center(2)) then
                    if (targ(3)>Current_box%Box_center(3)) then
                        call find_leaf_node(Current_box%Children(2)%BP,targ,p_box)
                    else
                        call find_leaf_node(Current_box%Children(6)%BP,targ,p_box)
                    endif
                else
                    if (targ(3)>Current_box%Box_center(3)) then
                        call find_leaf_node(Current_box%Children(4)%BP,targ,p_box)
                    else
                        call find_leaf_node(Current_box%Children(8)%BP,targ,p_box)
                    endif
                endif
            endif
        endif

return
end

subroutine eval_sigma_LRT(targ_vect,n_targ,FSS_1,sgma,Gx,Gy,Gz)
implicit none

    !List of calling arguments
    type ( Fast_Sigma_stuff ), pointer ::  FSS_1
    integer, intent(in) :: n_targ
    real ( kind = 8 ), intent(in) :: targ_vect(3,n_targ)
    real ( kind = 8 ), intent(out) :: sgma(n_targ),Gx(n_targ),Gy(n_targ),Gz(n_targ)

    !List of local variables
    real ( kind = 8 ) F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z,di2,my_exp,alpha,aux_center(3),sigma2_aux,lambda,beta
    integer count1,count2
    type ( TreeLRD ), pointer :: TreeLRD_1
    type ( Box ), pointer :: p_box



        TreeLRD_1 => FSS_1%TreeLRD_1
!        alpha=1.0d0/maxval(TreeLRD_1%W_sgmas_mem)**2/5.0d0
        alpha=FSS_1%alpha
        lambda=5.0d0
        beta=15.0d0
        do count1=1,n_targ
            call find_leaf_node(TreeLRD_1%Main_box,targ_vect(:,count1),p_box)
            F=0.0d0
            D=0.0d0
            dF_x=0.0d0
            dF_y=0.0d0
            dF_z=0.0d0
            dD_x=0.0d0
            dD_y=0.0d0
            dD_z=0.0d0
            do count2=1,56
                if (associated(p_box%Nearest(count2)%BP)) then
                    aux_center=p_box%Nearest(count2)%BP%Box_center
                    sigma2_aux=(p_box%Nearest(count2)%BP%Box_size/2)**10
                    di2=((aux_center(1)-targ_vect(1,count1))**2+(aux_center(2)&
                     &-targ_vect(2,count1))**2+(aux_center(3)-targ_vect(3,count1))**2)

                    my_exp=exp(-di2**5/(sigma2_aux))
                    F=F+p_box%Nearest(count2)%BP%Box_size/lambda*my_exp
                    !F=F+my_exp

                    dF_x=dF_x+p_box%Nearest(count2)%BP%Box_size/lambda*my_exp*(targ_vect(1,count1)-aux_center(1))/sigma2_aux
                    dF_y=dF_y+p_box%Nearest(count2)%BP%Box_size/lambda*my_exp*(targ_vect(2,count1)-aux_center(2))/sigma2_aux
                    dF_z=dF_z+p_box%Nearest(count2)%BP%Box_size/lambda*my_exp*(targ_vect(3,count1)-aux_center(3))/sigma2_aux
                    dD_x=dD_x+my_exp*(targ_vect(1,count1)-aux_center(1))/sigma2_aux
                    dD_y=dD_y+my_exp*(targ_vect(2,count1)-aux_center(2))/sigma2_aux
                    dD_z=dD_z+my_exp*(targ_vect(3,count1)-aux_center(3))/sigma2_aux
                    D=D+my_exp
                endif
                aux_center=p_box%Box_center
                sigma2_aux=(p_box%Box_size/beta)**2
                di2=((aux_center(1)-targ_vect(1,count1))**2+(aux_center(2)&
                &-targ_vect(2,count1))**2+(aux_center(3)-targ_vect(3,count1))**2)

                my_exp=exp(-di2/(2*sigma2_aux))
                F=F+p_box%Box_size/lambda*my_exp
                !F=F+my_exp
!                write (*,*) 'sigma target: ', p_box%Box_size/lambda

                dF_x=dF_x+p_box%Box_size/lambda*my_exp*(targ_vect(1,count1)-aux_center(1))/sigma2_aux
                dF_y=dF_y+p_box%Box_size/lambda*my_exp*(targ_vect(2,count1)-aux_center(2))/sigma2_aux
                dF_z=dF_z+p_box%Box_size/lambda*my_exp*(targ_vect(3,count1)-aux_center(3))/sigma2_aux
                dD_x=dD_x+my_exp*(targ_vect(1,count1)-aux_center(1))/sigma2_aux
                dD_y=dD_y+my_exp*(targ_vect(2,count1)-aux_center(2))/sigma2_aux
                dD_z=dD_z+my_exp*(targ_vect(3,count1)-aux_center(3))/sigma2_aux
                D=D+my_exp
            enddo
            Gx(count1)=(dF_x*D-F*dD_x)/D**2
            Gy(count1)=(dF_y*D-F*dD_y)/D**2
            Gz(count1)=(dF_z*D-F*dD_z)/D**2
            sgma(count1)=F/D
        enddo

return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!    END DEFINITION OF SUBROUTINES    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine s_interp_grad(x,y,yp)
implicit none

    real ( kind = 8 ), intent(in) :: x
    real ( kind = 8 ), intent(out) :: y, yp

        y=(erf((x-.5d0)/0.084d0)+1.0d0)/2.0d0
        yp=6.716542661282813d0*exp(-((x-0.5d0)/0.084d0)**2)

return
end



function s_interp(x)
implicit none

    real ( kind = 8 ), intent(in) :: x
    real ( kind = 8 ) :: s_interp

        s_interp=(erf((x-.5d0)/0.084d0)+1.0d0)/2.0d0

end function s_interp



function interp_3D(x_ent,y_ent,z_ent,P_vect) result(RR)
implicit none

    real ( kind = 8 ) RR
    real ( kind = 8 ), intent(in) :: P_vect(3,3,3)
    real ( kind = 8 ), intent(in) :: x_ent,y_ent,z_ent
    real ( kind = 8 ) c1pp,c3pp,cp1p,cp3p,cpp1,cpp3
    real ( kind = 8 ) a11p,a13p,a33p,a31p
    real ( kind = 8 ) a1p1,a1p3,a3p3,a3p1
    real ( kind = 8 ) ap11,ap13,ap33,ap31
!    real ( kind = 8 ) s_interp,interp_1D,interp_2D
    real ( kind = 8 ) x(2),y(2),z(2)

        x(1)=x_ent
        y(1)=y_ent
        z(1)=z_ent

        x(2)=s_interp(x(1))
        y(2)=s_interp(y(1))
        z(2)=s_interp(z(1))



        cpp1=interp_2D(x,y,P_vect(:,:,1))
        cpp3=interp_2D(x,y,P_vect(:,:,3))

        c1pp=interp_2D(y,z,P_vect(1,:,:))
        c3pp=interp_2D(y,z,P_vect(3,:,:))

        cp1p=interp_2D(x,z,P_vect(:,1,:))
        cp3p=interp_2D(x,z,P_vect(:,3,:))

        ap11=interp_1D(x,P_vect(:,1,1))
        ap13=interp_1D(x,P_vect(:,1,3))
        ap31=interp_1D(x,P_vect(:,3,1))
        ap33=interp_1D(x,P_vect(:,3,3))

        a1p1=interp_1D(y,P_vect(1,:,1))
        a1p3=interp_1D(y,P_vect(1,:,3))
        a3p1=interp_1D(y,P_vect(3,:,1))
        a3p3=interp_1D(y,P_vect(3,:,3))

        a11p=interp_1D(z,P_vect(1,1,:))
        a13p=interp_1D(z,P_vect(1,3,:))
        a31p=interp_1D(z,P_vect(3,1,:))
        a33p=interp_1D(z,P_vect(3,3,:))

        RR=cpp1*(1-z(2))+cpp3*z(2)
        RR=RR+c1pp*(1-x(2))+c3pp*x(2)
        RR=RR+cp1p*(1-y(2))+cp3p*y(2)

        RR=RR-ap11*(1-y(2))*(1-z(2))
        RR=RR-ap13*(1-y(2))*z(2)
        RR=RR-ap31*y(2)*(1-z(2))
        RR=RR-ap33*y(2)*z(2)

        RR=RR-a1p1*(1-x(2))*(1-z(2))
        RR=RR-a1p3*(1-x(2))*z(2)
        RR=RR-a3p1*x(2)*(1-z(2))
        RR=RR-a3p3*x(2)*z(2)

        RR=RR-a11p*(1-x(2))*(1-y(2))
        RR=RR-a13p*(1-x(2))*y(2)
        RR=RR-a31p*x(2)*(1-y(2))
        RR=RR-a33p*x(2)*y(2)

        RR=RR+P_vect(1,1,1)*(1-x(2))*(1-y(2))*(1-z(2))
        RR=RR+P_vect(1,1,3)*(1-x(2))*(1-y(2))*z(2)
        RR=RR+P_vect(1,3,1)*(1-x(2))*y(2)*(1-z(2))
        RR=RR+P_vect(1,3,3)*(1-x(2))*y(2)*z(2)

        RR=RR+P_vect(3,1,1)*x(2)*(1-y(2))*(1-z(2))
        RR=RR+P_vect(3,1,3)*x(2)*(1-y(2))*z(2)
        RR=RR+P_vect(3,3,1)*x(2)*y(2)*(1-z(2))
        RR=RR+P_vect(3,3,3)*x(2)*y(2)*z(2)

end function interp_3D



function interp_2D(u,v,P_vect) result(RR)
implicit none

    real ( kind = 8 ) RR
    real ( kind = 8 ), intent(in) :: P_vect(3,3)
    real ( kind = 8 ), intent(in) :: u(2),v(2)

    real ( kind = 8 ) u_aux(2),v_aux(2),a,b,c,d
    real ( kind = 8 ) P_aux(2,2)
!    real ( kind = 8 ) s_interp,interp_1D


        if (P_vect(2,2)>0.0d0) then
            if (u(1)<0.5d0) then
                if (v(1)<0.5d0) then
                    u_aux(1)=2*u(1)
                    v_aux(1)=2*v(1)
                    P_aux=P_vect(1:2,1:2)
                else
                    u_aux(1)=2*u(1)
                    v_aux(1)=2*(v(1)-0.5d0)
                    P_aux=P_vect(1:2,2:3)
                endif
            else
                if (v(1)<0.5d0) then
                    u_aux(1)=2*(u(1)-0.5d0)
                    v_aux(1)=2*v(1)
                    P_aux=P_vect(2:3,1:2)
                else
                    u_aux(1)=2*(u(1)-0.5d0)
                    v_aux(1)=2*(v(1)-0.5d0)
                    P_aux=P_vect(2:3,2:3)
                endif
            endif
            u_aux(2)=s_interp(u_aux(1))
            v_aux(2)=s_interp(v_aux(1))

            RR=P_aux(1,1)*(1-u_aux(2))*(1-v_aux(2))
            RR=RR+P_aux(2,1)*u_aux(2)*(1-v_aux(2))
            RR=RR+P_aux(1,2)*(1-u_aux(2))*v_aux(2)
            RR=RR+P_aux(2,2)*u_aux(2)*v_aux(2)
        else
            a=interp_1D(u,P_vect(:,1))
            b=interp_1D(u,P_vect(:,3))
            c=interp_1D(v,P_vect(1,:))
            d=interp_1D(v,P_vect(3,:))
            RR=a*(1-v(2))+b*v(2)
            RR=RR+c*(1-u(2))+d*u(2)
            RR=RR-P_vect(1,1)*(1-u(2))*(1-v(2))
            RR=RR-P_vect(3,1)*u(2)*(1-v(2))
            RR=RR-P_vect(1,3)*(1-u(2))*v(2)
            RR=RR-P_vect(3,3)*u(2)*v(2)
        endif

end function interp_2D


function interp_1D(t,P_vect) result(RR)
implicit none

    real ( kind = 8 ) RR
    real ( kind = 8 ), intent(in) :: P_vect(3)
    real ( kind = 8 ), intent(in) :: t(2)
!    real ( kind = 8 ) s_interp
    real ( kind = 8 ) t_aux(2)

        if (P_vect(2)<0.0d0) then
            RR=P_vect(1)*(1-t(2))+P_vect(3)*t(2)
        else
            if (t(1)<0.5) then
                t_aux(1)=2*t(1)
                t_aux(2)=s_interp(t_aux(1))
                RR=P_vect(1)*(1-t_aux(2))+P_vect(2)*t_aux(2)
            else
                t_aux(1)=2*(t(1)-0.5d0)
                t_aux(2)=s_interp(t_aux(1))

                RR=P_vect(2)*(1-t_aux(2))+P_vect(3)*t_aux(2)
            endif
        endif

end function interp_1D






subroutine interp_1D_grad(t,P_vect,RR,RRp)
implicit none

    real ( kind = 8 ), intent(out) :: RR,RRp
    real ( kind = 8 ), intent(in) :: P_vect(3)
    real ( kind = 8 ), intent(in) :: t(3)
!    real ( kind = 8 ) s_interp
    real ( kind = 8 ) t_aux(3)

        if (P_vect(2)<0.0d0) then
            RR=P_vect(1)*(1-t(2))+P_vect(3)*t(2)
            RRp=P_vect(1)*(-t(3))+P_vect(3)*t(3)
        else
            if (t(1)<0.5) then
                t_aux(1)=2*t(1)
                call s_interp_grad(t_aux(1),t_aux(2),t_aux(3))
                RR=P_vect(1)*(1-t_aux(2))+P_vect(2)*t_aux(2)
                RRp=P_vect(1)*(-t_aux(3))+P_vect(2)*t_aux(3)
            else
                t_aux(1)=2*(t(1)-0.5d0)
                call s_interp_grad(t_aux(1),t_aux(2),t_aux(3))
                RR=P_vect(2)*(1-t_aux(2))+P_vect(3)*t_aux(2)
                RRp=P_vect(2)*(-t_aux(3))+P_vect(3)*t_aux(3)
            endif
        endif

return
end

subroutine interp_2D_grad(u,v,P_vect,RR,RRu,RRv)
implicit none

    real ( kind = 8 ), intent(out) :: RR,RRu,RRv
    real ( kind = 8 ), intent(in) :: P_vect(3,3)
    real ( kind = 8 ), intent(in) :: u(3),v(3)

    real ( kind = 8 ) u_aux(3),v_aux(3)
    real ( kind = 8 ) a,b,c,d,au,av,bu,bv,cu,cv,du,dv
    real ( kind = 8 ) P_aux(2,2)
!    real ( kind = 8 ) s_interp,interp_1D


        if (P_vect(2,2)>0.0d0) then
            if (u(1)<0.5d0) then
                if (v(1)<0.5d0) then
                    u_aux(1)=2*u(1)
                    v_aux(1)=2*v(1)
                    P_aux=P_vect(1:2,1:2)
                else
                    u_aux(1)=2*u(1)
                    v_aux(1)=2*(v(1)-0.5d0)
                    P_aux=P_vect(1:2,2:3)
                endif
            else
                if (v(1)<0.5d0) then
                    u_aux(1)=2*(u(1)-0.5d0)
                    v_aux(1)=2*v(1)
                    P_aux=P_vect(2:3,1:2)
                else
                    u_aux(1)=2*(u(1)-0.5d0)
                    v_aux(1)=2*(v(1)-0.5d0)
                    P_aux=P_vect(2:3,2:3)
                endif
            endif
            call s_interp_grad(u_aux(1),u_aux(2),u_aux(3))
            call s_interp_grad(v_aux(1),v_aux(2),v_aux(3))

            RR=P_aux(1,1)*(1-u_aux(2))*(1-v_aux(2))
            RR=RR+P_aux(2,1)*u_aux(2)*(1-v_aux(2))
            RR=RR+P_aux(1,2)*(1-u_aux(2))*v_aux(2)
            RR=RR+P_aux(2,2)*u_aux(2)*v_aux(2)

            RRu=P_aux(1,1)*(-u_aux(3))*(1-v_aux(2))
            RRu=RRu+P_aux(2,1)*u_aux(3)*(1-v_aux(2))
            RRu=RRu+P_aux(1,2)*(-u_aux(3))*v_aux(2)
            RRu=RRu+P_aux(2,2)*u_aux(3)*v_aux(2)

            RRv=P_aux(1,1)*(1-u_aux(2))*(-v_aux(3))
            RRv=RRv+P_aux(2,1)*u_aux(2)*(-v_aux(3))
            RRv=RRv+P_aux(1,2)*(1-u_aux(2))*v_aux(3)
            RRv=RRv+P_aux(2,2)*u_aux(2)*v_aux(3)

        else

            call interp_1D_grad(u,P_vect(:,1),a,au)
            call interp_1D_grad(u,P_vect(:,3),b,bu)
            call interp_1D_grad(v,P_vect(1,:),c,cv)
            call interp_1D_grad(v,P_vect(3,:),d,dv)

            RR=a*(1-v(2))+b*v(2)
            RR=RR+c*(1-u(2))+d*u(2)
            RR=RR-P_vect(1,1)*(1-u(2))*(1-v(2))
            RR=RR-P_vect(3,1)*u(2)*(1-v(2))
            RR=RR-P_vect(1,3)*(1-u(2))*v(2)
            RR=RR-P_vect(3,3)*u(2)*v(2)

            RRu=au*(1-v(2))+bu*v(2)
            RRu=RRu+c*(-u(3))+d*u(3)
            RRu=RRu-P_vect(1,1)*(-u(3))*(1-v(2))
            RRu=RRu-P_vect(3,1)*u(3)*(1-v(2))
            RRu=RRu-P_vect(1,3)*(-u(3))*v(2)
            RRu=RRu-P_vect(3,3)*u(3)*v(2)

            RRv=a*(-v(3))+b*v(3)
            RRv=RRv+cv*(1-u(2))+dv*u(2)
            RRv=RRv-P_vect(1,1)*(1-u(2))*(-v(3))
            RRv=RRv-P_vect(3,1)*u(2)*(-v(3))
            RRv=RRv-P_vect(1,3)*(1-u(2))*v(3)
            RRv=RRv-P_vect(3,3)*u(2)*v(3)

        endif

return
end

subroutine interp_3D_grad(x_ent,y_ent,z_ent,P_vect,RR,RRx,RRy,RRz)
implicit none

    real ( kind = 8 ),intent(out) :: RR,RRx,RRy,RRz
    real ( kind = 8 ), intent(in) :: P_vect(3,3,3)
    real ( kind = 8 ), intent(in) :: x_ent,y_ent,z_ent
    real ( kind = 8 ) c1pp,c3pp,cp1p,cp3p,cpp1,cpp3
    real ( kind = 8 ) c1ppy,c1ppz,c3ppy,c3ppz,cp1px,cp1pz
    real ( kind = 8 ) cp3px,cp3pz,cpp1x,cpp1y,cpp3x,cpp3y
    real ( kind = 8 ) a11p,a11pz,a13p,a13pz,a33p,a33pz,a31p,a31pz
    real ( kind = 8 ) a1p1,a1p3,a3p3,a3p1
    real ( kind = 8 ) a1p1y,a1p3y,a3p3y,a3p1y
    real ( kind = 8 ) ap11,ap13,ap33,ap31
    real ( kind = 8 ) ap11x,ap13x,ap33x,ap31x
!    real ( kind = 8 ) s_interp,interp_1D,interp_2D
    real ( kind = 8 ) x(3),y(3),z(3)

        x(1)=x_ent
        y(1)=y_ent
        z(1)=z_ent

        call s_interp_grad(x(1),x(2),x(3))
        call s_interp_grad(y(1),y(2),y(3))
        call s_interp_grad(z(1),z(2),z(3))

        call interp_2D_grad(x,y,P_vect(:,:,1),cpp1,cpp1x,cpp1y)
        call interp_2D_grad(x,y,P_vect(:,:,3),cpp3,cpp3x,cpp3y)

        call interp_2D_grad(y,z,P_vect(1,:,:),c1pp,c1ppy,c1ppz)
        call interp_2D_grad(y,z,P_vect(3,:,:),c3pp,c3ppy,c3ppz)

        call interp_2D_grad(x,z,P_vect(:,1,:),cp1p,cp1px,cp1pz)
        call interp_2D_grad(x,z,P_vect(:,3,:),cp3p,cp3px,cp3pz)

        call interp_1D_grad(x,P_vect(:,1,1),ap11,ap11x)
        call interp_1D_grad(x,P_vect(:,1,3),ap13,ap13x)
        call interp_1D_grad(x,P_vect(:,3,1),ap31,ap31x)
        call interp_1D_grad(x,P_vect(:,3,3),ap33,ap33x)

        call interp_1D_grad(y,P_vect(1,:,1),a1p1,a1p1y)
        call interp_1D_grad(y,P_vect(1,:,3),a1p3,a1p3y)
        call interp_1D_grad(y,P_vect(3,:,1),a3p1,a3p1y)
        call interp_1D_grad(y,P_vect(3,:,3),a3p3,a3p3y)

        call interp_1D_grad(z,P_vect(1,1,:),a11p,a11pz)
        call interp_1D_grad(z,P_vect(1,3,:),a13p,a13pz)
        call interp_1D_grad(z,P_vect(3,1,:),a31p,a31pz)
        call interp_1D_grad(z,P_vect(3,3,:),a33p,a33pz)

        RR=cpp1*(1-z(2))+cpp3*z(2)
        RR=RR+c1pp*(1-x(2))+c3pp*x(2)
        RR=RR+cp1p*(1-y(2))+cp3p*y(2)

        RR=RR-ap11*(1-y(2))*(1-z(2))
        RR=RR-ap13*(1-y(2))*z(2)
        RR=RR-ap31*y(2)*(1-z(2))
        RR=RR-ap33*y(2)*z(2)

        RR=RR-a1p1*(1-x(2))*(1-z(2))
        RR=RR-a1p3*(1-x(2))*z(2)
        RR=RR-a3p1*x(2)*(1-z(2))
        RR=RR-a3p3*x(2)*z(2)

        RR=RR-a11p*(1-x(2))*(1-y(2))
        RR=RR-a13p*(1-x(2))*y(2)
        RR=RR-a31p*x(2)*(1-y(2))
        RR=RR-a33p*x(2)*y(2)

        RR=RR+P_vect(1,1,1)*(1-x(2))*(1-y(2))*(1-z(2))
        RR=RR+P_vect(1,1,3)*(1-x(2))*(1-y(2))*z(2)
        RR=RR+P_vect(1,3,1)*(1-x(2))*y(2)*(1-z(2))
        RR=RR+P_vect(1,3,3)*(1-x(2))*y(2)*z(2)

        RR=RR+P_vect(3,1,1)*x(2)*(1-y(2))*(1-z(2))
        RR=RR+P_vect(3,1,3)*x(2)*(1-y(2))*z(2)
        RR=RR+P_vect(3,3,1)*x(2)*y(2)*(1-z(2))
        RR=RR+P_vect(3,3,3)*x(2)*y(2)*z(2)



        RRx=cpp1x*(1-z(2))+cpp3x*z(2)
        RRx=RRx+c1pp*(-x(3))+c3pp*x(3)
        RRx=RRx+cp1px*(1-y(2))+cp3px*y(2)

        RRx=RRx-ap11x*(1-y(2))*(1-z(2))
        RRx=RRx-ap13x*(1-y(2))*z(2)
        RRx=RRx-ap31x*y(2)*(1-z(2))
        RRx=RRx-ap33x*y(2)*z(2)

        RRx=RRx-a1p1*(-x(3))*(1-z(2))
        RRx= RRx-a1p3*(-x(3))*z(2)
        RRx= RRx-a3p1*x(3)*(1-z(2))
        RRx= RRx-a3p3*x(3)*z(2)

        RRx= RRx-a11p*(-x(3))*(1-y(2))
        RRx= RRx-a13p*(-x(3))*y(2)
        RRx= RRx-a31p*x(3)*(1-y(2))
        RRx= RRx-a33p*x(3)*y(2)

        RRx= RRx+P_vect(1,1,1)*(-x(3))*(1-y(2))*(1-z(2))
        RRx= RRx+P_vect(1,1,3)*(-x(3))*(1-y(2))*z(2)
        RRx= RRx+P_vect(1,3,1)*(-x(3))*y(2)*(1-z(2))
        RRx= RRx+P_vect(1,3,3)*(-x(3))*y(2)*z(2)

        RRx= RRx+P_vect(3,1,1)*x(3)*(1-y(2))*(1-z(2))
        RRx= RRx+P_vect(3,1,3)*x(3)*(1-y(2))*z(2)
        RRx= RRx+P_vect(3,3,1)*x(3)*y(2)*(1-z(2))
        RRx= RRx+P_vect(3,3,3)*x(3)*y(2)*z(2)


        RRy=cpp1y*(1-z(2))+cpp3y*z(2)
        RRy=RRy+c1ppy*(1-x(2))+c3ppy*x(2)
        RRy=RRy+cp1p*(-y(3))+cp3p*y(3)

        RRy=RRy-ap11*(-y(3))*(1-z(2))
        RRy=RRy-ap13*(-y(3))*z(2)
        RRy=RRy-ap31*y(3)*(1-z(2))
        RRy=RRy-ap33*y(3)*z(2)

        RRy=RRy-a1p1y*(1-x(2))*(1-z(2))
        RRy=RRy-a1p3y*(1-x(2))*z(2)
        RRy=RRy-a3p1y*x(2)*(1-z(2))
        RRy=RRy-a3p3y*x(2)*z(2)

        RRy=RRy-a11p*(1-x(2))*(-y(3))
        RRy=RRy-a13p*(1-x(2))*y(3)
        RRy=RRy-a31p*x(2)*(-y(3))
        RRy=RRy-a33p*x(2)*y(3)

        RRy=RRy+P_vect(1,1,1)*(1-x(2))*(-y(3))*(1-z(2))
        RRy=RRy+P_vect(1,1,3)*(1-x(2))*(-y(3))*z(2)
        RRy=RRy+P_vect(1,3,1)*(1-x(2))*y(3)*(1-z(2))
        RRy=RRy+P_vect(1,3,3)*(1-x(2))*y(3)*z(2)

        RRy=RRy+P_vect(3,1,1)*x(2)*(-y(3))*(1-z(2))
        RRy=RRy+P_vect(3,1,3)*x(2)*(-y(3))*z(2)
        RRy=RRy+P_vect(3,3,1)*x(2)*y(3)*(1-z(2))
        RRy=RRy+P_vect(3,3,3)*x(2)*y(3)*z(2)


        RRz=cpp1*(-z(3))+cpp3*z(3)
        RRz=RRz+c1ppz*(1-x(2))+c3ppz*x(2)
        RRz=RRz+cp1pz*(1-y(2))+cp3pz*y(2)

        RRz=RRz-ap11*(1-y(2))*(-z(3))
        RRz=RRz-ap13*(1-y(2))*z(3)
        RRz=RRz-ap31*y(2)*(-z(3))
        RRz=RRz-ap33*y(2)*z(3)

        RRz=RRz-a1p1*(1-x(2))*(-z(3))
        RRz=RRz-a1p3*(1-x(2))*z(3)
        RRz=RRz-a3p1*x(2)*(-z(3))
        RRz=RRz-a3p3*x(2)*z(3)

        RRz=RRz-a11pz*(1-x(2))*(1-y(2))
        RRz=RRz-a13pz*(1-x(2))*y(2)
        RRz=RRz-a31pz*x(2)*(1-y(2))
        RRz=RRz-a33pz*x(2)*y(2)

        RRz=RRz+P_vect(1,1,1)*(1-x(2))*(1-y(2))*(-z(3))
        RRz=RRz+P_vect(1,1,3)*(1-x(2))*(1-y(2))*z(3)
        RRz=RRz+P_vect(1,3,1)*(1-x(2))*y(2)*(-z(3))
        RRz=RRz+P_vect(1,3,3)*(1-x(2))*y(2)*z(3)

        RRz=RRz+P_vect(3,1,1)*x(2)*(1-y(2))*(-z(3))
        RRz=RRz+P_vect(3,1,3)*x(2)*(1-y(2))*z(3)
        RRz=RRz+P_vect(3,3,1)*x(2)*y(2)*(-z(3))
        RRz=RRz+P_vect(3,3,3)*x(2)*y(2)*z(3)


return
end







end module Mod_Fast_Sigma
