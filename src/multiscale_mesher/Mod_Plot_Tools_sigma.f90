Module Mod_Plot_Tools_sigma
  use Mod_TreeLRD
  use ModType_Smooth_Surface
  use Mod_Fast_Sigma
  implicit none

public :: plot_sigma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!    DEFINITION OF SUBROUTINES     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

subroutine plot_sigma(FSS_1, Geometry1,adapt_flag)
implicit none

!List of calling arguments
type ( Fast_Sigma_stuff ), pointer :: FSS_1
type ( Geometry ), pointer :: Geometry1
integer *8, intent(in) :: adapt_flag

!List of local variables
type ( TreeLRD ), pointer :: TreeLRD_1
character ( len=100 ) nombre,filename,plot_name
real ( kind = 8 ) x_min,x_max,y_min,y_max,z_min,z_max,Lx,Ly,Lz
real ( kind = 8 ) xp_min,xp_max,yp_min,yp_max,zp_min,zp_max
real ( kind = 8 ), allocatable :: F_plot(:,:),targ_vect(:,:),sgma(:),sgma_x(:),sgma_y(:),sgma_z(:)
integer *8 N_plot,M_plot,count,count1,count2,icount,n_targ


    TreeLRD_1 => FSS_1%TreeLRD_1

!    do count=1,Geometry1%ntri_sk
!        write (*,*) TreeLRD_1%W_Pts_mem(1,count),TreeLRD_1%W_Pts_mem(2,count),TreeLRD_1%W_Pts_mem(3,count)
!        write (*,*) TreeLRD_1%W_sgmas_mem(count)
!    enddo
    plot_name='./plot_tools/plot_points'

call plot_curve_3D(TreeLRD_1%W_Pts_mem(1,:),TreeLRD_1%W_Pts_mem(2,:),TreeLRD_1%W_Pts_mem(3,:),TreeLRD_1%ntot_W_Pts_mem,plot_name)

    write(*,*) 'inside2'

    M_plot=800
    N_plot=800

    allocate(F_plot(M_plot,N_plot))
    allocate(targ_vect(3,M_plot*N_plot))
    allocate(sgma(M_plot*N_plot))
    allocate(sgma_x(M_plot*N_plot))
    allocate(sgma_y(M_plot*N_plot))
    allocate(sgma_z(M_plot*N_plot))

    x_min=minval(TreeLRD_1%W_Pts_mem(1,:))
    x_max=maxval(TreeLRD_1%W_Pts_mem(1,:))
    y_min=minval(TreeLRD_1%W_Pts_mem(2,:))
    y_max=maxval(TreeLRD_1%W_Pts_mem(2,:))
    z_min=minval(TreeLRD_1%W_Pts_mem(3,:))
    z_max=maxval(TreeLRD_1%W_Pts_mem(3,:))
    Lx=x_max-x_min
    Ly=y_max-y_min
    Lz=z_max-z_min

    xp_max=x_max+Lx/4.0d0
    xp_min=x_min-Lx/4.0d0
    yp_max=y_max+Ly/4.0d0
    yp_min=y_min-Ly/4.0d0
    zp_max=z_max+Lz/4.0d0
    zp_min=z_min-Lz/4.0d0


    icount=1
    do count2=1,M_plot
        do count1=1,N_plot
!!! Plane ZY
            targ_vect(3,icount)=zp_min+(zp_max-zp_min)*(count1-1.0d0)/M_plot
            targ_vect(2,icount)=yp_min+(yp_max-yp_min)*(count2-1.0d0)/N_plot
            targ_vect(1,icount)=(x_max+x_min)/2.0d0+(x_max-x_min)/4.0d0*0.0d0
!            call initial_guess_sgma(TreeLRD_1%Main_box,targ_vect(:,icount),sgma(icount))
            icount=icount+1
        enddo
    enddo
    n_targ=N_plot*M_plot
    call function_eval_sigma(FSS_1,targ_vect,n_targ,sgma,sgma_x,sgma_y,sgma_z,adapt_flag)
    write (*,*) 'sgma_val: ',sgma(1)
    plot_name='./plot_tools/plot_sigma_1'
    call plot2D_v3(targ_vect(1,:),targ_vect(2,:),targ_vect(3,:),sgma,N_plot,M_plot,plot_name)


    icount=1
    do count2=1,M_plot
        do count1=1,N_plot
!!! Plane XY
            targ_vect(1,icount)=xp_min+(xp_max-xp_min)*(count1-1.0d0)/M_plot
            targ_vect(2,icount)=yp_min+(yp_max-yp_min)*(count2-1.0d0)/N_plot
            targ_vect(3,icount)=(z_max+z_min)/2.0d0+(z_max-z_min)/4.0d0*0.0d0+.50d0
!            call initial_guess_sgma(TreeLRD_1%Main_box,targ_vect(:,icount),sgma(icount))
            icount=icount+1
        enddo
    enddo
    n_targ=N_plot*M_plot
    call function_eval_sigma(FSS_1,targ_vect,n_targ,sgma,sgma_x,sgma_y,sgma_z,adapt_flag)

    plot_name='./plot_tools/plot_sigma_2'
    call plot2D_v3(targ_vect(1,:),targ_vect(2,:),targ_vect(3,:),sgma,N_plot,M_plot,plot_name)

    icount=1
    do count2=1,M_plot
        do count1=1,N_plot
!!! Plane XZ
            targ_vect(1,icount)=xp_min+(xp_max-xp_min)*(count1-1.0d0)/M_plot
            targ_vect(3,icount)=zp_min+(zp_max-zp_min)*(count2-1.0d0)/N_plot
            targ_vect(2,icount)=(y_max+y_min)/2.0d0
!            call initial_guess_sgma(TreeLRD_1%Main_box,targ_vect(:,icount),sgma(icount))
            icount=icount+1
        enddo
    enddo
    n_targ=N_plot*M_plot

    call function_eval_sigma(FSS_1,targ_vect,n_targ,sgma,sgma_x,sgma_y,sgma_z,adapt_flag)

    plot_name='./plot_tools/plot_sigma_3'
    call plot2D_v3(targ_vect(1,:),targ_vect(2,:),targ_vect(3,:),sgma,N_plot,M_plot,plot_name)

    deallocate(F_plot)
    deallocate(targ_vect)
    deallocate(sgma)
    deallocate(sgma_x)
    deallocate(sgma_y)
    deallocate(sgma_z)

    call plot_tree_tool(TreeLRD_1)

return
end


subroutine plot_tree_tool(TreeLRD_1)
implicit none

!List of calling arguments
type ( TreeLRD ), pointer :: TreeLRD_1

!List of local variables
character ( len=100 ) plot_name
real ( kind = 8 ), allocatable :: W_boxes(:)
integer *8 icount,tot_num_box

    plot_name='./plot_tools/plot_tree'
    tot_num_box=total_number_leaf_boxes(TreeLRD_1%Main_box)
    allocate(W_boxes(tot_num_box*4))
    icount=1
    call pile_leaf_boxes(TreeLRD_1%Main_box,W_boxes,icount,tot_num_box)
    call plot_tree(W_boxes,tot_num_box,TreeLRD_1%W_Pts_mem,TreeLRD_1%ntot_W_Pts_mem,plot_name)
    deallocate(W_boxes)

return
end


recursive subroutine pile_boxes(Current_box,W_boxes,icount,tot_num_box)
implicit none

!List of calling arguments
type ( Box ), pointer :: Current_box
integer *8, intent(in) :: tot_num_box!List of local variables
integer *8, intent(inout) :: icount!List of local variables
real ( kind = 8 ), intent(inout) :: W_boxes(4*tot_num_box)

!List of local variables
integer *8 count1

    if (Current_box%is_leaf) then
        W_boxes(icount)=Current_box%Box_center(1)
        W_boxes(icount+1)=Current_box%Box_center(2)
        W_boxes(icount+2)=Current_box%Box_center(3)
        W_boxes(icount+3)=Current_box%Box_size
        icount=icount+4
    else
        do count1=1,8
            call pile_boxes(Current_box%Children(count1)%BP,W_boxes,icount,tot_num_box)
        enddo
        W_boxes(icount)=Current_box%Box_center(1)
        W_boxes(icount+1)=Current_box%Box_center(2)
        W_boxes(icount+2)=Current_box%Box_center(3)
        W_boxes(icount+3)=Current_box%Box_size
        icount=icount+4
    endif

return
end


recursive subroutine pile_leaf_boxes(Current_box,W_boxes,icount,tot_num_box)
implicit none

!List of calling arguments
type ( Box ), pointer :: Current_box
integer *8, intent(in) :: tot_num_box!List of local variables
integer *8, intent(inout) :: icount!List of local variables
real ( kind = 8 ), intent(inout) :: W_boxes(4*tot_num_box)

!List of local variables
integer *8 count1

    if (Current_box%is_leaf) then
        W_boxes(icount)=Current_box%Box_center(1)
        W_boxes(icount+1)=Current_box%Box_center(2)
        W_boxes(icount+2)=Current_box%Box_center(3)
        W_boxes(icount+3)=Current_box%Box_size
        icount=icount+4
    else
        do count1=1,8
            call pile_leaf_boxes(Current_box%Children(count1)%BP,W_boxes,icount,tot_num_box)
        enddo
    endif

return
end


end module Mod_Plot_Tools_sigma
