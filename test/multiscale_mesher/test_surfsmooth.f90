
program smoother

  implicit none

  integer :: N, count,nrefine, ifplot
  integer :: adapt_flag, ifflatten
  integer :: interp_flag
  integer :: norder_skel, norder_smooth

  character (len=100) :: fnamein, fnameout_root,name_aux
  character (len=21) :: plot_name
  character (len=8) :: istr1,istr2
  character (len=2) :: arg_comm

  real *8 rlam
  integer i, ier, ifiletype
  


  call prini(6,13)


  ! order with which to discretize the skeleton patches (pick
  ! something high-order)
  norder_skel = 12

  ! order with which to discretize the smooth patches, choose
  ! something reasonable: 4, 6, 8, 10, etc.
  norder_smooth = 4

  ! Define number of refinements of smooth surface to be output in
  ! the go3 format
  ! 
  !
  nrefine = 2
  ! nrefine=1

  ! this is to enable adaptativity (otherwise sigma is constant)
  ! adapt_flag = 0  ->  no adaptivity, mean triangle size
  ! adapt_flag = 1  ->  some adaptivity, alpha form
  ! adapt_flag = 2  ->  full recursive definition, slightly slower
  adapt_flag = 1


  !
  !
  ! rlam flag decides the proportionality value for \sigma_{j}
  ! in relation to triangle diameter
  ! \sigma_{j} = D_{j}/rlam
  !
  rlam = 10.0d0 !(usual value)

  !rlam = .5d0
  !rlam = 1
  !rlam = 2.5d0


!
! specify the msh file to read in
!

!    fnamein='../../geometries/meshes/cuboid_a1_b2_c1p3.tri'
!    fnameout_root='../../geometries/cuboid_a1_b2_c1p3'

!    fnamein='../../geometries/meshes/prism_50.gidmsh'
!    fnameout_root='../../geometries/prism_50'

    fnamein='../../geometries/meshes/sphere.msh'
    fnameout_root='../../geometries/sphere'
    
!    fnamein = '../../geometries/meshes/cow_new.msh'
!    fnameout_root = '../../geometries/cow_new'

!    fnamein = '../../geometries/meshes/cow_new_gmshv4.msh'
!    fnameout_root = '../../geometries/cow_new'

!    fnamein = '../../geometries/meshes/lens_r00.msh'
!    fnameout_root = '../../geometries/lens_r00'

!    fnamein = '../../geometries/meshes/lens_r00_gmshv4.msh'
!    fnameout_root = '../../geometries/lens_r00'
!


!    Determine file type
!        * ifiletype = 1, for .msh
!        * ifiletype = 2, for .tri
!        * ifiletype = 3, for .gidmsh
!        * ifiletype = 4, for .msh gmsh v2
!        * ifiletype = 5, for .msh gmsh v4
    ier = 0
    call get_filetype(fnamein, ifiletype, ier)

    if(ier.ne.0) then
      print *, "File type not recognized"
      stop
    endif

    ier = 0
    call multiscale_mesher_unif_refine(fnamein, ifiletype, norder_skel, &
       norder_smooth, nrefine, adapt_flag, rlam, fnameout_root, ier)

end program smoother

