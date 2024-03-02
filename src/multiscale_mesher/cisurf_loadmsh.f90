
!
! (c) 2019 Felipe Vico and Michael O'Neil
! oneil@cims.nyu.edu
!
! This collection of subroutines loads a msh file, puts info into the
! Geometry variable, and records the resulting smooth surface in a
! *.gov file
!

subroutine readgeometry(Geometry1, filename, norder_skel, &
    norder_smooth)
  use ModType_Smooth_Surface
  implicit none
  integer ifiletype, ier
  type (Geometry) :: Geometry1
  character(len=*) :: filename
  integer :: n_order_sf, norder_skel, norder_smooth

  !
  ! This subroutine open a msh file and load the information in a
  ! variable of type Geometry
  !
  ! Input
  !   filename - name of file
  !   n_order_skel - discretization order for representing surface
  !   norder_smooth - discretization order for computing integrals
  !

  if (norder_skel .gt. 20) then
    call prinf('norder_skel too large = *', norder_skel, 1)
    stop
  end if

  if (norder_smooth .gt. 20) then
    call prinf('norder_smooth too large = *', norder_smooth, 1)
    stop
  end if

!    Determine file type
!        * ifiletype = 1, for .msh
!        * ifiletype = 2, for .tri
!        * ifiletype = 3, for .gidmsh
!        * ifiletype = 4, for .gmsh
  call get_filetype(filename, ifiletype, ier)

  if(ier.ne.0) then
    print *, "File type not recognized"
    stop
  endif


  if (ifiletype.eq.1) then
    call readmsh(Geometry1, filename, norder_skel, norder_smooth)
  elseif (ifiletype.eq.2) then
    call readtri(Geometry1, filename, norder_skel, norder_smooth)
  elseif (ifiletype.eq.3) then
    call readgidmsh(Geometry1, filename, norder_skel, norder_smooth)
  elseif (ifiletype.eq.4) then
    call read_gmsh_v2(Geometry1, filename, norder_skel, norder_smooth)
  elseif (ifiletype.eq.5) then
     print *, "mesh type not supported"
     stop
!    call read_gmsh_v2(Geometry1, filename, norder_skel, norder_smooth)
  else
    write (*,*) 'Mesh type not supported' 
    stop
  endif

  return
end subroutine readgeometry
!
!
!
!
!



subroutine readmsh(Geometry1, filename, norder_skel, norder_smooth)
  use ModType_Smooth_Surface
  implicit none

  !! This subroutine open a msh file and load the information in a
  !! variable of type Geometry

  !
  ! Input
  !   filename - the file to read
  !   norder_skel - order to discretize the skeleton patches
  !   norder_smooth - order to discretize the smoothed patches
  !
  ! Output
  !   Geometry1 - data structure for geometry
  !

  !List of calling arguments
  type (Geometry), intent(inout) :: Geometry1     !! where the geometry will be loaded
  character(len=*), intent(in) :: filename         !! name of the msh file
  integer :: n_order_sf, nsk, nsf
  integer  :: norder_skel, norder_smooth
  
  integer umio,i,m,N,j,aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
  integer :: ierror


  Geometry1%ifflat = 0

  open(UNIT=8, FILE=trim(filename), STATUS='OLD', ACTION='READ', IOSTAT=ierror)
  read(8,*) aux1,aux2,aux3,m, N

  
  nsk = (norder_skel+1)*(norder_skel+2)/2
  
  Geometry1%norder_skel = norder_skel
  Geometry1%nskel = nsk

  !Geometry1%norder_smooth = norder_smooth
  !Geometry1%nsmooth = n_order_sf

  Geometry1%norder_smooth = norder_smooth
  nsf = (norder_smooth+1)*(norder_smooth+2)/2
  Geometry1%nsmooth = nsf

  call prinf('num points on smooth mesh = *', nsf*n, 1)
  
  Geometry1%n_order_sf = nsf
  Geometry1%npoints=m
  Geometry1%ntri=N
  Geometry1%ntri_sk=N
  Geometry1%n_Sf_points=N*nsf
  Geometry1%n_Sk_points=N*nsk

  if (allocated(Geometry1%Points)) then
    deallocate(Geometry1%Points)
  endif
  if (allocated(Geometry1%Tri)) then
    deallocate(Geometry1%Tri)
  endif

  allocate(Geometry1%Points(3,m))
  allocate(Geometry1%Tri(6,N))

  do j=1,m
    read(8,*) aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
    read(8,*) Geometry1%Points(1,j),Geometry1%Points(2,j),Geometry1%Points(3,j)
  enddo

  read(8,*) aux1

  do j=1,N
    read(8,*) aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
    read(8,*) aux1,aux2,aux3,aux4,aux5,aux6
    Geometry1%Tri(1,j)=aux1
    Geometry1%Tri(2,j)=aux2
    Geometry1%Tri(3,j)=aux3
    Geometry1%Tri(4,j)=aux4
    Geometry1%Tri(5,j)=aux5
    Geometry1%Tri(6,j)=aux6
  enddo
  
  close (8)

  return
end subroutine readmsh





subroutine readgidmsh(Geometry1, filename, norder_skel, norder_smooth)
  use ModType_Smooth_Surface
  implicit none

  !! This subroutine open a msh file and load the information in a
  !! variable of type Geometry

  !
  ! Input
  !   filename - the file to read
  !   norder_skel - order to discretize the skeleton patches
  !   norder_smooth - order to discretize the smoothed patches
  !
  ! Output
  !   Geometry1 - data structure for geometry
  !

  !List of calling arguments
  type (Geometry), intent(inout) :: Geometry1     !! where the geometry will be loaded
  character(len=*), intent(in) :: filename         !! name of the msh file
  character(len=1000) :: tmp1
  integer :: n_order_sf, nsk, nsf
  integer  :: norder_skel, norder_smooth
  
  integer :: umio,i,m,N,j,aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
  integer :: node, nnodes, maxnodes
  integer, allocatable :: elems(:,:)
  integer :: ielem, nelems, maxelems
  double precision :: x, y, z, d, dmin
  double precision, allocatable :: xs(:), ys(:), zs(:)
  integer :: ierror, iflag


  Geometry1%ifflat = 0

  open(UNIT=8, FILE=trim(filename), STATUS='OLD', ACTION='READ', &
     IOSTAT=ierror)


  read(8,'(a)') tmp1
  read(8,'(a)') tmp1

  ! count the number of nodes and elements
  !
  nnodes = 0
  iflag = 0
  do while (iflag.eq.0) 
    read (8,'(a)') tmp1
    if (index(tmp1, 'End Coordinates') > 0) then
      iflag = 1
    else
      nnodes = nnodes + 1
    endif
  enddo

  do i = 1,100
    read(8,*) tmp1
    if (index(tmp1, 'Elements') > 0) exit
  end do

  iflag = 0
  nelems = 0

  do while (iflag.eq.0)
    read (8,'(a)') tmp1
    if (index(tmp1, 'End Elements') > 0) then
      iflag = 1
    else
      nelems = nelems + 1 
    endif
  enddo
  rewind(8)

  n = nelems

  Geometry1%norder_skel = norder_skel
  nsk = (norder_skel+1)*(norder_skel+2)/2
  Geometry1%nskel = nsk

  Geometry1%norder_smooth = norder_smooth
  nsf = (norder_smooth+1)*(norder_smooth+2)/2
  Geometry1%nsmooth = nsf

  Geometry1%n_order_sf = nsf

  m = nnodes
  Geometry1%npoints=m
  Geometry1%ntri=n
  Geometry1%ntri_sk=n
  Geometry1%n_Sf_points=n*nsf
  Geometry1%n_Sk_points=n*nsk

  if (allocated(Geometry1%Points)) then
    deallocate(Geometry1%Points)
  endif
  if (allocated(Geometry1%Tri)) then
    deallocate(Geometry1%Tri)
  endif

  allocate(Geometry1%Points(3,m))
  allocate(Geometry1%Tri(6,n))

  read(8,'(a)') tmp1
  read(8,'(a)') tmp1

  do i=1,nnodes
    read(8,*) node, Geometry1%Points(1,i), Geometry1%Points(2,i), &
      Geometry1%Points(3,i)
  enddo

  
  do i = 1,100
    read(8,*) tmp1
    if (index(tmp1, 'Elements') > 0) exit
  end do

  do i=1,nelems
    read (8, *) ielem, Geometry1%Tri(1,i), Geometry1%Tri(2,i), &
    Geometry1%Tri(3,i), Geometry1%Tri(4,i), &
    Geometry1%Tri(5,i), Geometry1%Tri(6,i)
  enddo

  close(8)
  

  return
end subroutine readgidmsh
!
!
!
!
!

subroutine readtri(Geometry1,filename, norder_skel, norder_smooth)
  use ModType_Smooth_Surface
  implicit none

  !! This subroutine open a msh file and load the information in a
  !! variable of type Geometry

  type (Geometry), intent(inout) :: Geometry1     !! where the geometry will be loaded
  character(len=*), intent(in) :: filename         !! name of the msh file
  integer :: norder_smooth, norder_skel

  integer umio,i,m,N,j,aux1,aux2,aux3,ipointer
  integer :: ierror, nsk,nsf

  ! set the flag for flat vs quadratic
  Geometry1%ifflat = 1

  
  open(UNIT=8, FILE=trim(filename), STATUS='OLD', ACTION='READ', &
    IOSTAT=ierror)
  read(8,*) m, N


  call prinf('npoints = *', m, 1)
  call prinf('ntri = *', n, 1)


  nsk = (norder_skel+1)*(norder_skel+2)/2
  call prinf('nsk = *', nsk, 1)
  !stop

  Geometry1%norder_skel = norder_skel
  Geometry1%nskel = nsk
  
  Geometry1%norder_smooth = norder_smooth
  nsf = (norder_smooth+1)*(norder_smooth+2)/2
  
  Geometry1%n_order_sf=nsf
  Geometry1%nsmooth = nsf
  Geometry1%npoints=m+N*3
  Geometry1%ntri=N
  Geometry1%ntri_sk=N
  Geometry1%n_Sf_points=N*nsf
  Geometry1%n_Sk_points=N*nsk

  if (allocated(Geometry1%Points)) then
    deallocate(Geometry1%Points)
  endif
  if (allocated(Geometry1%Tri)) then
    deallocate(Geometry1%Tri)
  endif
  allocate(Geometry1%Points(3,Geometry1%npoints))
  allocate(Geometry1%Tri(6,N))

  do j=1,m
    read(8,*) Geometry1%Points(1,j),Geometry1%Points(2,j), &
       Geometry1%Points(3,j)
  enddo


  ipointer=m+1

  do j=1,N
    read(8,*) aux1,aux2,aux3
    Geometry1%Tri(1,j)=aux1
    Geometry1%Tri(2,j)=aux2
    Geometry1%Tri(3,j)=aux3
    Geometry1%Tri(4,j)=ipointer
    Geometry1%Tri(5,j)=ipointer+1
    Geometry1%Tri(6,j)=ipointer+2
    Geometry1%Points(:,ipointer) = & 
      (Geometry1%Points(:,Geometry1%Tri(1,j)) + &
       Geometry1%Points(:,Geometry1%Tri(2,j)))/2.0d0
    Geometry1%Points(:,ipointer+1) = &
      (Geometry1%Points(:,Geometry1%Tri(2,j)) + &
       Geometry1%Points(:,Geometry1%Tri(3,j)))/2.0d0
    Geometry1%Points(:,ipointer+2) = &
      (Geometry1%Points(:,Geometry1%Tri(1,j)) + &
       Geometry1%Points(:,Geometry1%Tri(3,j)))/2.0d0
    ipointer=ipointer+3
  enddo
  close (8)

  return
end subroutine readtri




subroutine read_gmsh_v2(Geometry1, filename, norder_skel, norder_smooth)
  use ModType_Smooth_Surface
  implicit none

  ! This subroutine opens a v2 gmsh file and load the information in a
  ! variable of type Geometry. Currently only reads second order 
  ! quads and triangles. 
  !
  ! Todo: Support for first order stuff 

  !
  ! Input
  !   filename - the file to read
  !   norder_skel - order to discretize the skeleton patches
  !   norder_smooth - order to discretize the smoothed patches
  !
  ! Output
  !   Geometry1 - data structure for geometry
  !

  !List of calling arguments
  type (Geometry), intent(inout) :: Geometry1     !! where the geometry will be loaded
  character(len=*), intent(in) :: filename         !! name of the msh file
  character(len=1000) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7
  character(len=1000) :: cline
  integer :: n_order_sf, nsk, nsf
  integer  :: norder_skel, norder_smooth
  
  integer :: umio,i,m,N,j,aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
  integer :: node, nnodes, maxnodes
  integer, allocatable :: elements(:,:), element(:)
  integer :: ielem, nelems, maxelems
  double precision :: x, y, z, d, dmin
  double precision, allocatable :: xs(:), ys(:), zs(:)
  integer :: ierror,iunit,korder,kpols,itype
  integer :: io,numnodes,ind,numelem,nel,ntri,ntag,lll
  integer :: itype_tri, itype_quad, kpols_quad, kpols_tri


  Geometry1%ifflat = 0

  iunit = 899

  open(UNIT=iunit, FILE=trim(filename), STATUS='OLD', ACTION='READ', IOSTAT=ierror)


  korder = 2
  kpols_tri = 6
  kpols_quad = 9
  kpols = (korder+1)*(korder+2)/2

  itype_tri = 9
  itype_quad = 10

!
!   Read the number of points
!
!

  do
    read(iunit, *, iostat=io) cline

    if (io .ne. 0) exit

    if (trim(cline) .eq. '$Nodes') then
      print *, 'Reading nodes . . . '
      read(iunit,*) numnodes
      print *, 'Number of nodes = ', numnodes
      print *
      
      allocate(xs(numnodes),ys(numnodes),zs(numnodes))
      do i = 1,numnodes
        read (iunit,*) ind, x, y, z
        xs(i) = x
        ys(i) = y
        zs(i) = z
      end do
      exit

    end if
  enddo

  nel = 6 
  ntri = 0

  do
    read(iunit, *, iostat=io) cline
    if (io .ne. 0) exit

    if (trim(cline) .eq. '$Elements') then
      print *, 'Reading elements . . . '
      read(iunit,*) numelem
      print *, 'Number of elements = ', numelem
      print *


      ntri = 0
      do i = 1,numelem
        
        read(iunit, '(a)', iostat=io) cline
        read (cline,*) ind, ielem, ntag

        if (ielem .eq. itype_tri) ntri = ntri + 1
        if (ielem .eq. itype_quad) ntri = ntri + 2
      enddo
      exit
    endif
   enddo


  rewind(899)
  n = ntri

  Geometry1%norder_skel = norder_skel
  nsk = (norder_skel+1)*(norder_skel+2)/2
  Geometry1%nskel = nsk


  Geometry1%norder_smooth = norder_smooth
  nsf = (norder_smooth+1)*(norder_smooth+2)/2
  Geometry1%nsmooth = nsf

  Geometry1%n_order_sf = nsf

  m = numnodes
  Geometry1%npoints=m
  Geometry1%ntri=N
  Geometry1%ntri_sk=N
  Geometry1%n_Sf_points=N*nsf
  Geometry1%n_Sk_points=N*nsk

  if (allocated(Geometry1%Points)) then
    deallocate(Geometry1%Points)
  endif
  if (allocated(Geometry1%Tri)) then
    deallocate(Geometry1%Tri)
  endif

  allocate(Geometry1%Points(3,m))
  allocate(Geometry1%Tri(6,N))

  do j=1,m
    Geometry1%Points(1,j) = xs(j)
    Geometry1%Points(2,j) = ys(j)
    Geometry1%Points(3,j) = zs(j)
  enddo

   ntri = 0
   do
    read(iunit, *, iostat=io) cline
    if (io .ne. 0) exit

    if (trim(cline) .eq. '$Elements') then
      read(iunit,*) numelem
      do i = 1, numelem

        read(iunit, '(a)', iostat=io) cline
        read (cline,*) ind, ielem, ntag
        if(ielem.eq.itype_tri) lll= 1+1+1+ntag+nel
        if(ielem.eq.itype_quad) lll= 14 
        allocate(element(lll))
        read(cline,*) element
        if (ielem.eq.itype_tri) then
          ntri = ntri + 1
          do j = 1,6
            Geometry1%Tri(j,ntri) = element(j+3+ntag)
          end do

        end if

        if (ielem.eq.itype_quad) then
          ntri = ntri + 1
          Geometry1%Tri(1,ntri) = element(6)
          Geometry1%Tri(2,ntri) = element(7)
          Geometry1%Tri(3,ntri) = element(9)
          Geometry1%Tri(4,ntri) = element(10)
          Geometry1%Tri(5,ntri) = element(14)
          Geometry1%Tri(6,ntri) = element(13)


          ntri = ntri + 1
          Geometry1%Tri(1,ntri) = element(8)
          Geometry1%Tri(2,ntri) = element(9)
          Geometry1%Tri(3,ntri) = element(7)
          Geometry1%Tri(4,ntri) = element(12)
          Geometry1%Tri(5,ntri) = element(14)
          Geometry1%Tri(6,ntri) = element(11)

        endif
        deallocate(element)
        
      end do

    end if
    
  end do

  close(iunit)

  return
end subroutine read_gmsh_v2
!
!
!
!
!

subroutine record_Geometry(Geometry1,filename)
  use ModType_Smooth_Surface
  implicit none

  !List of calling arguments
  type (Geometry), intent(in) :: Geometry1
  character (len=*) filename

  !List of local variables
  integer umio,count1,count2,flag,norder_smooth
  integer :: ierror

  open(8, FILE=trim(filename),STATUS='REPLACE')
  norder_smooth = Geometry1%norder_smooth
  
  write(8,*) norder_smooth
  write(8,*) Geometry1%ntri
!  write(8,*) Geometry1%n_Sf_points
  do count1=1,Geometry1%n_Sf_points
    write(8,*) Geometry1%S_smooth(1,count1)
  enddo
  do count1=1,Geometry1%n_Sf_points
    write(8,*) Geometry1%S_smooth(2,count1)
  enddo
  do count1=1,Geometry1%n_Sf_points
    write(8,*) Geometry1%S_smooth(3,count1)
  enddo

  do count1=1,Geometry1%n_Sf_points
	write(8,*) Geometry1%du_smooth(1,count1)
  enddo
  do count1=1,Geometry1%n_Sf_points
	write(8,*) Geometry1%du_smooth(2,count1)
  enddo
  do count1=1,Geometry1%n_Sf_points
	write(8,*) Geometry1%du_smooth(3,count1)
  enddo
  
  do count1=1,Geometry1%n_Sf_points
	write(8,*) Geometry1%dv_smooth(1,count1)
  enddo
  do count1=1,Geometry1%n_Sf_points
	write(8,*) Geometry1%dv_smooth(2,count1)
  enddo
  do count1=1,Geometry1%n_Sf_points
	write(8,*) Geometry1%dv_smooth(3,count1)
  enddo

  do count1=1,Geometry1%n_Sf_points
    write(8,*) Geometry1%N_smooth(1,count1)
  enddo
  do count1=1,Geometry1%n_Sf_points
    write(8,*) Geometry1%N_smooth(2,count1)
  enddo
  do count1=1,Geometry1%n_Sf_points
    write(8,*) Geometry1%N_smooth(3,count1)
  enddo

  close (8)

  return
end subroutine record_Geometry



subroutine get_filetype(filename, ifiletype, ier)
!
!  This subroutine determines the mesh type of a given
!  file. Currently supported formats include
!
!  * .msh
!  * .tri
!  * .gidmsh
!  * .msh (from gmsh v2 or v4)
!
!  Input arguments:
!    - filename: string
!         file name
!  Output arguments
!    - ifiletype: integer
!        * ifiletype = 1, for .msh
!        * ifiletype = 2, for .tri
!        * ifiletype = 3, for .gidmsh
!        * ifiletype = 4, for .gmsh
!    - ier: integer
!        Error code
!        * ier = 0, successful execution
!        * ier = 8, file format not recognized
!
!
  implicit real *8 (a-h,o-z)
  integer i1,i2, io
  character (len=*) filename
  character (len=100) fstr
  character *100, fstr_gid_tritest, fstr_gmshtest

  ier = 0
  fstr_gid_tritest = 'MESH dimension 3 ElemType Triangle Nnode '
  fstr_gmshtest = '$MeshFormat'

! 
  open(unit=33, file=trim(filename), status='old')

!  Check if it is a .msh file
  io = 0   
  read(33,*, iostat = io) i1,i2, i3, i4, i5
  if(io.eq.0) then
    ifiletype = 1
    return
  endif
      
!   Check if it is a .tri file now      
  rewind(33)
  read(33,*, iostat = io) i1,i2
  if(io.eq.0) then
    ifiletype = 2
    return
  endif

!
!  Now check if it is a gidmsh
!
  rewind(33)
  read(33,'(a)', iostat = io) fstr
      
  len1 = len(trim(fstr_gid_tritest))
  len2 = len(trim(fstr))
  if(trim(fstr(1:len1)) .eq. trim(fstr_gid_tritest)) then
     read(fstr((len1+1):len2),*,iostat=io) iind
     ifiletype = 3
     return
  endif

  len1 = len(trim(fstr_gmshtest))
  if(trim(fstr(1:len1)) .eq. trim(fstr_gmshtest)) then
     read(33,*,iostat=io) tmp, tmp2, tmp3
     if(abs(tmp-2).le.1) then
       ifiletype = 4
     elseif(abs(tmp-4).le.1) then
       ifiletype = 5
     else
       print *, "Invalid gmsh file type"
       ier = 8
     endif
     return
  endif

  print *, "Invalid file type"
  ier = 8
  close(33)

  return
  end
