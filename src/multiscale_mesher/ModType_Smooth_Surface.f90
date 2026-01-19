!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!                    MODULE MAIN INFORMATION                          !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! DESCRIPTION:
!!!
!!! This module contains types necessary to smooth a quadratic surface
!!! The main type introduced is Geometry, that contains all relevant information about
!!! the skeleton, the Newton search pseudonormals and store the resulting smooth surface
!!!
!!! INDEX OF SUBROUTINES AND FUNCTIONS: (none)
!!!
!!! INDEX OF DERIVED TYPES:
!!!
!!! Public:
!!!
!!! 1-Geometry              !! Main type. Contains all the information required for runing newton
!!! 2-Variable_Matrix       !! simple variable lenght matrix (used to compute the 'normal on a vertex')
!!! 3-My_cell               !! simple group of Variable_Matrix (used to compute the 'normal on a vertex')
!!!
!!! Private: (none)
!!!
!!! MODULES REQUIRED: (none)
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module ModType_Smooth_Surface
implicit none

!! This type contains all the information necessary to smooth a surface starting
!! from a quadratic patch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!    DEFINITION OF TYPES     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type, public :: Geometry
    real ( kind = 8 ), allocatable :: S_smooth(:,:)             !Points on the real smooth surface
    real ( kind = 8 ), allocatable :: N_smooth(:,:)             !Normals on the real smooth surface
    real ( kind = 8 ), allocatable :: skeleton_Points(:,:)      !Integration nodes on the skeleton to compute F (F=0.5 is the surface)
    real ( kind = 8 ), allocatable :: skeleton_w(:)             !Integration weights on the skeleton to compute F (F=0.5 is the surface)
    real ( kind = 8 ), allocatable :: skeleton_N(:,:)           !Normal Vector on the skeleton at the nodes (to compute the double layer)
    real ( kind = 8 ), allocatable :: ru_smooth(:,:)            !u vectro on the real smooth surface
    real ( kind = 8 ), allocatable :: rv_smooth(:,:)            !v vectro on the real smooth surface
	real ( kind = 8 ), allocatable :: du_smooth(:,:)            !u vectro on the real smooth surface
    real ( kind = 8 ), allocatable :: dv_smooth(:,:)            !v vectro on the real smooth surface
    real ( kind = 8 ), allocatable :: ds_smooth(:)            !u vectro on the real smooth surface
    real ( kind = 8 ), allocatable :: w_smooth(:)               !Integration weigths on the real smooth surface
!

    ! That was removed from here. There will be only one copy of the centroids and sgmas at centroids in the tree.
!    real ( kind = 8 ), allocatable :: Centroids(:,:)            !Centroids of each triangle of the skeleton (to compute sgma(x))
!    real ( kind = 8 ), allocatable :: sgmas(:)                  !Values of sgma on each centroid, proportional to the side of the triangle
    !

    ! flag for triangle type, flat or quadratic
    integer *8 :: ifflat
    
    !Triangles of the skeleton (each triangle 6 points)
    integer *8, allocatable :: Tri(:,:)

    !Points that define the msh file (Each quadratic triangle has 6 points)
    real ( kind = 8 ), allocatable :: Points(:,:)

    !Pseudo-normals defined on each vertex of each triangle of the skeleton
    real ( kind = 8 ), allocatable :: Normal_Vert(:,:)          
    real ( kind = 8 ), allocatable :: Base_Points(:,:)          !base base points of eachs mooth point on the skeleton
    real ( kind = 8 ), allocatable :: Base_Points_N(:,:)        !Pseudo-normals defined on each base point of the surface
    real ( kind = 8 ), allocatable :: Base_Points_U(:,:)        !U vector defined on each base point of the smooth surface
    real ( kind = 8 ), allocatable :: Base_Points_V(:,:)        !U vector defined on each base point of the smooth surface
    real ( kind = 8 ), allocatable :: Dummy_targ(:,:)        !Dummy targets on pseudonormals
    real ( kind = 8 ), allocatable :: height(:)                 !Value of the root in newton where F(h)=1/2 for each point on the smooth surface
    integer *8 n_dummy_targ   !total number of dummy targets

!Edges and Boundary is use to compute the line integral in the non-FMM case (beyond that line integral, the kernel is 1/r and erf=1)
    integer *8, allocatable :: Edges(:,:)             !Set of 3 integers with the location of the points Points(1:3,:) for each edge in the skeleton (set of quadratic triangles)
    real ( kind = 8 ), allocatable :: skelet_line_p(:,:,:)            ! integration nodes along each boundary of each triangle on the skeleton
    real ( kind = 8 ), allocatable :: skelet_line_dl(:,:,:)           ! vector dl along each boundary to compute line integrals (3,n_nodes per side,edge number)
    integer *8, allocatable :: Boundary(:,:)          ! Set of 3 Edges for each triangle in the skeleton:  Boundary(n_edge=1 to 3,num of triangle)
!!
    integer *8 npoints                                !Total number of points in the skeleton
    integer *8 n_Sf_points                            !total number of points on the real smooth surface
    integer *8 n_Sk_points                            !total number of integration nodes on the skeleton
    integer *8 ntri                                   !Total number of triangles on the smooth surface
    integer *8 ntri_sk                                !Total number of
    !triangles on the skeleton

    ! order of smooth discretization, and points per triangle
    integer *8 norder_skel
    integer *8 nskel
    
    integer *8 norder_smooth
    integer *8 nsmooth
    !number of nodes per smooth triangle (45 or 78)
    integer *8 n_order_sf   

end type Geometry

type, public :: Variable_Matrix
    integer *8 n_Mat, current_n_Mat
    real ( kind = 8 ), allocatable :: Mat(:,:)
end type Variable_Matrix

type, public :: My_cell
    type (Variable_Matrix), allocatable :: Var_Mat(:)
    integer *8 n_Cell
end type My_cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!    END DEFINITION OF TYPES     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ModType_Smooth_Surface
