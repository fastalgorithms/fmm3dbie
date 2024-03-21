subroutine multiscale_mesher_unif_refine(fnamein, ifiletype, norder_skel, &
   norder_smooth, nrefine, adapt_sigma, rlam, fnameout_root, ier)
!
!  Given an input flat/second order triangulated mesh 
!  specified in .gidmsh, .msh, gmsh v2, or .tri formats,
!  run the mutiscale surface smoother to obtain a high
!  order discretization of the surface
!
!  NOTE: in order to determine the file type of your 
!  input file, you can use the utility
!  get_filetype, see cisurf_loadmsh.f90 for documentation
!  
!  Input arguments:
!    * fnamein: string
!        Input File name
!    * ifiletype: integer
!        * ifiletype = 1, for .msh
!        * ifiletype = 2, for .tri
!        * ifiletype = 3, for .gidmsh
!        * ifiletype = 4, for .msh gmsh v2
!        * ifiletype = 5, for .msh gmsh v4
!    * norder_skel: integer
!        Order of discretization on skeleton mesh for 
!        computing integrals for the level-set
!    * norder_smooth: integer
!        order of discretization of the surface
!        on output
!    * nrefine: integer
!        number of refinements of the skeleton mesh
!    * adapt_sigma: integer
!        flag for choosing uniform or adaptive mollifier
!        adapt_sigma = 0, use uniform mollifier
!        adapt_sigma = 1, use multiscale mollifier
!    * rlam: real *8
!        rlam decides the proportionality value for \sigma_{j}
!        in relation to triangle diameter 
!        \sigma_{j} = D_{j}/rlam. The \sigma_{j} is then used
!        as a variance for the Gaussian centered at triangle
!        centroid for obtaining target dependent variance of
!        the final mollifier.
!
!        Small values of rlam lead to less smoothing, larger
!        values of rlam lead to more smoothing. Recommended
!        range of parameters: 2.5 - 10
!    * fnameout_root: string
!        root name for storing output files
!
!  Output arguments:
!    * ier: error code
!         ier = 0, implies successful execution
!         ier = 2, implies error reading input mesh file
!         ier = 3, implies newton failed to converge for one of
!                  the refinements
!
!
   

  use ModType_Smooth_Surface
  use Mod_TreeLRD
  use Mod_Plot_Tools_sigma
  use Mod_Fast_Sigma
  use Mod_Feval
  use Mod_Smooth_Surface

  implicit none

  type ( Geometry ), pointer :: Geometry1 => null ()
  type ( Feval_stuff ), pointer :: Feval_stuff_1 => null ()

  integer :: N, count,nrefine, ifplot
  integer :: adapt_sigma, ifflatten
  integer :: interp_flag
  integer :: norder_skel, norder_smooth
  integer :: ifiletype

  character (len=*) :: fnamein, fnameout_root

  character(:), allocatable :: fname_aux
  character (len=21) :: plot_name
  character (len=8) :: istr1,istr2
  character (len=2) :: arg_comm
  double precision :: x0,y0,z0, ra
  double precision, allocatable :: time_report(:), error_report(:)
  double precision :: err_skel,rlam,t1,t2,omp_get_wtime

  integer i, ll, len1

  integer ier

  

  allocate(Geometry1)
  allocate(error_report(nrefine+1))


  !
  ! specify the msh file to read in
  !

  ier = 0
  ! load in the msh file
  call readgeometry(Geometry1, trim(fnamein), ifiletype, norder_skel, &
      norder_smooth, ier)
  if(ier.ne.0) then
    print *, "Error reading input mesh file"
    return
  endif

  print *, "successfully read input mesh"
  print *, "-----------------------------"
  print *, ""
  print *, "Starting surface smoother now"
  print *
  count = 0
  print *, "-----------------------------------"
  print *, "Refinement:", count
  print *, "-----------------------------------"


  ! dump out discretization points on the skeleton mesh
  call funcion_skeleton(Geometry1)
  call funcion_normal_vert(Geometry1)


!
!  Compute centroid
!
  x0 = 0
  y0 = 0
  z0 = 0
  ra = 0
  do i=1,Geometry1%n_Sk_points
    x0 = x0 + Geometry1%skeleton_Points(1,i)*Geometry1%skeleton_w(i)
    y0 = y0 + Geometry1%skeleton_Points(2,i)*Geometry1%skeleton_w(i)
    z0 = z0 + Geometry1%skeleton_Points(3,i)*Geometry1%skeleton_w(i)
    ra = ra + Geometry1%skeleton_w(i)
  enddo
  x0 = x0/ra
  y0 = y0/ra
  z0 = z0/ra


  call start_Feval_tree(Feval_stuff_1, Geometry1, rlam)
  call funcion_Base_Points(Geometry1)
  call find_smooth_surface(Geometry1, Feval_stuff_1, adapt_sigma, ier)

  if(ier.ne.0) then
    print *, "Newton failed to converge at refinement:",  count
    return
  endif

  len1 = len(trim(fnameout_root))

  ll = len1 + 12
  allocate(character(ll) :: fname_aux)

  count = 0

  write(istr1,"(I2.2)") count
  write(istr2,"(I2.2)") norder_smooth
  fname_aux = trim(fnameout_root)//'_o'//trim(istr2)// &
     '_r'//trim(istr1)//'.go3'
  call record_Geometry(Geometry1,fname_aux)

!
!  check Gauss error
!
  call check_Gauss(Geometry1,x0,y0,z0,error_report(1))
  print *, " "
  print *, " "
  print *, "-----------------------------------"
  do count=1,nrefine
    print *, "Refinement:", count 
    print *, "-----------------------------------"
     
    call refine_geometry_smart(Geometry1)

    call funcion_Base_Points(Geometry1)

    call find_smooth_surface(Geometry1,Feval_stuff_1,adapt_sigma, ier)
    if(ier.ne.0) then
      print *, "Newton failed to converge at refinement:",  count
      return
    endif

    write(istr1,"(I2.2)") count
    write(istr2,"(I2.2)") norder_smooth
    fname_aux = trim(fnameout_root)//'_o'//trim(istr2)// &
       '_r'//trim(istr1)//'.go3'
    call record_Geometry(Geometry1,fname_aux)

    call check_Gauss(Geometry1,x0,y0,z0,error_report(count+1))
    print *, " "
    print *, " "
    print *, "-----------------------------------"
  enddo
  
  print *, ""
  print *, ""
  print *, "Report for Error in Gauss' theorem"
  print *, "-----------------------------------"
  do count=0,nrefine
    write (*,*) 'Refinement nº: ',int(count,4), '  Error: ', &
        real(error_report(count+1),4)
  enddo

end subroutine multiscale_mesher_unif_refine 
!
!
!
!
!
!
subroutine multiscale_mesher_unif_refine_cfname(fnamein, ifiletype, norder_skel, &
   norder_smooth, nrefine, adapt_sigma, rlam, fnameout_root, ier)
!
!  Same routine as the above multiscale mesher, the only difference
!  is that the filename on input are of type cstring which are converted
!  to fortran strings and then sent to the fortran version of the routine
!
!  Given an input flat/second order triangulated mesh 
!  specified in .gidmsh, .msh, gmsh v2, or .tri formats,
!  run the mutiscale surface smoother to obtain a high
!  order discretization of the surface
!
!  NOTE: in order to determine the file type of your 
!  input file, you can use the utility
!  get_filetype, see cisurf_loadmsh.f90 for documentation
!  
!  Input arguments:
!    * fnamein: cstring 
!        Input File name
!    * ifiletype: integer
!        * ifiletype = 1, for .msh
!        * ifiletype = 2, for .tri
!        * ifiletype = 3, for .gidmsh
!        * ifiletype = 4, for .msh gmsh v2
!        * ifiletype = 5, for .msh gmsh v4
!    * norder_skel: integer
!        Order of discretization on skeleton mesh for 
!        computing integrals for the level-set
!    * norder_smooth: integer
!        order of discretization of the surface
!        on output
!    * nrefine: integer
!        number of refinements of the skeleton mesh
!    * adapt_sigma: integer
!        flag for choosing uniform or adaptive mollifier
!        adapt_sigma = 0, use uniform mollifier
!        adapt_sigma = 1, use multiscale mollifier
!    * rlam: real *8
!        rlam decides the proportionality value for \sigma_{j}
!        in relation to triangle diameter 
!        \sigma_{j} = D_{j}/rlam. The \sigma_{j} is then used
!        as a variance for the Gaussian centered at triangle
!        centroid for obtaining target dependent variance of
!        the final mollifier.
!
!        Small values of rlam lead to less smoothing, larger
!        values of rlam lead to more smoothing. Recommended
!        range of parameters: 2.5 - 10
!    * fnameout_root: cstring
!        root name for storing output files
!
!  Output arguments:
!    * ier: error code
!         ier = 0, implies successful execution
!         ier = 2, implies error reading input mesh file
!         ier = 3, implies newton failed to converge for one of
!                  the refinements
!
!
   

  use iso_c_binding

  implicit none


  integer :: N, count,nrefine, ifplot
  integer :: adapt_sigma, ifflatten
  integer :: interp_flag
  integer :: norder_skel, norder_smooth
  integer :: ifiletype

  character (kind=c_char), dimension(*) :: fnamein, fnameout_root
  character (len=:), allocatable :: fortran_fnamein, fortran_fnameout
  real *8 :: rlam
  integer :: ier


  integer ilen, i

  ilen = 0
  do 
    if (fnamein(ilen+1) == C_NULL_CHAR) exit
    ilen = ilen + 1
  enddo

  allocate(character(len=ilen) :: fortran_fnamein)
  fortran_fnamein = transfer(fnamein(1:ilen),fortran_fnamein)

  ilen = 0
  do 
    if (fnameout_root(ilen+1) == C_NULL_CHAR) exit
    ilen = ilen + 1
  enddo

  allocate(character(len=ilen) :: fortran_fnameout)
  fortran_fnameout = transfer(fnameout_root(1:ilen),fortran_fnameout)

  call multiscale_mesher_unif_refine(fortran_fnamein, ifiletype, norder_skel, &
   norder_smooth, nrefine, adapt_sigma, rlam, fortran_fnameout, ier)

  
  return
end subroutine multiscale_mesher_unif_refine_cfname 
!
!
!
!
!
!
subroutine get_filetype_cfname(fnamein, ifiletype, ier)
  use iso_c_binding  
  implicit none
  character (kind=c_char), dimension(*) :: fnamein 
  character (len=:), allocatable :: fortran_fnamein
  integer ilen, i
  integer ifiletype, ier

  ilen = 0
  do 
    if (fnamein(ilen+1) == C_NULL_CHAR) exit
    ilen = ilen + 1
  enddo

  allocate(character(len=ilen) :: fortran_fnamein)
  fortran_fnamein = transfer(fnamein(1:ilen),fortran_fnamein)

  call get_filetype(fortran_fnamein, ifiletype, ier)


  return
end subroutine get_filetype_cfname
