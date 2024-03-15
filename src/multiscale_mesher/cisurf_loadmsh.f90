
!
! (c) 2019 Felipe Vico and Michael O'Neil
! oneil@cims.nyu.edu
!
! This collection of subroutines loads a msh file, puts info into the
! Geometry variable, and records the resulting smooth surface in a
! *.gov file
!

subroutine readgeometry(Geometry1, filename, ifiletype, norder_skel, &
    norder_smooth, ier)
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
!  Input arguments:
!    - filename: string
!         file name
!    - ifiletype: integer
!        * ifiletype = 1, for .msh
!        * ifiletype = 2, for .tri
!        * ifiletype = 3, for .gidmsh
!        * ifiletype = 4, for .msh gmsh v2
!        * ifiletype = 5, for .msh gmsh v4
!   n_order_skel - discretization order for representing surface
!   norder_smooth - discretization order for computing integrals
!
! Output
!   ier: error code
!     ier = 0, successful execution
!     ier = 2, trouble in reading file
!     ier = 4, norder = norder_skel too large
!     ier = 5, nquad = norder_smooth too large
!

  ier = 0
  if (norder_skel .gt. 20) then
    call prinf('norder_skel too large = *', norder_skel, 1)
    ier = 4
    return
  end if

  if (norder_smooth .gt. 20) then
    call prinf('norder_smooth too large = *', norder_smooth, 1)
    ier = 5
    return
  end if



  if (ifiletype.eq.1) then
    call readmsh(Geometry1, filename, norder_skel, norder_smooth, ier)
  elseif (ifiletype.eq.2) then
    call readtri(Geometry1, filename, norder_skel, norder_smooth, ier)
  elseif (ifiletype.eq.3) then
    call readgidmsh(Geometry1, filename, norder_skel, norder_smooth, ier)
  elseif (ifiletype.eq.4) then
    call read_gmsh_v2(Geometry1, filename, norder_skel, norder_smooth, ier)
  elseif (ifiletype.eq.5) then
    call read_gmsh_v4(Geometry1, filename, norder_skel, norder_smooth, ier)
  else
    write (*,*) 'Mesh type not supported'
    ier = 1
    return
  endif

  return
end subroutine readgeometry
!
!
!
!
!



subroutine readmsh(Geometry1, filename, norder_skel, norder_smooth, ier)
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
  !   Geometry1 - data structure for geometry (if successful execution)
  !   ier - error code, ier = 0, implies successful execution
  !

  !List of calling arguments
  type (Geometry), intent(inout) :: Geometry1     !! where the geometry will be loaded
  character(len=*), intent(in) :: filename         !! name of the msh file
  integer :: n_order_sf, nsk, nsf
  integer  :: norder_skel, norder_smooth
  
  integer umio,i,m,N,j,aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
  integer :: ierror, ierror1, ierror2
  integer :: ier


  Geometry1%ifflat = 0
 
  ier = 0
  open(UNIT=8, FILE=trim(filename), STATUS='OLD', ACTION='READ', IOSTAT=ierror)
  if(ierror.ne.0) then
    ier = 2
    return
  endif
  read(8,*, iostat=ierror) aux1,aux2,aux3,m, N
  if(ierror.ne.0) then
    ier = 2
    return
  endif
  
  
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

  ierror = 0
  do j=1,m
    read(8,*, iostat=ierror1) aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
    read(8,*, iostat=ierror2) Geometry1%Points(1,j), &
       Geometry1%Points(2,j),Geometry1%Points(3,j)
    ierror = ierror1 + ierror2
  enddo

  read(8,*, iostat=ierror1) aux1
  ierror = ierror1 + ierror2

  do j=1,N
    read(8,*, iostat=ierror1) aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
    read(8,*, iostat=ierror2) aux1,aux2,aux3,aux4,aux5,aux6
    Geometry1%Tri(1,j)=aux1
    Geometry1%Tri(2,j)=aux2
    Geometry1%Tri(3,j)=aux3
    Geometry1%Tri(4,j)=aux4
    Geometry1%Tri(5,j)=aux5
    Geometry1%Tri(6,j)=aux6

    ierror = ierror1 + ierror2
  enddo
  if(ierror.ne.0) then
    ier = 2
    return
  endif
  
  close (8)

  return
end subroutine readmsh





subroutine readgidmsh(Geometry1, filename, norder_skel, norder_smooth, ier)
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
  integer :: ierror1, ierror2
  integer :: ier


  Geometry1%ifflat = 0

  ierror = 0
  ier = 0
  open(UNIT=8, FILE=trim(filename), STATUS='OLD', ACTION='READ', &
     IOSTAT=ierror)
  if(ierror.ne.0) then
    ier = 2
    return
  endif


  read(8,'(a)', iostat=ierror1) tmp1
  read(8,'(a)', iostat=ierror2) tmp1

  if(ierror1+ierror2.gt.0) then
    ier = 2
    return

  endif

  ! count the number of nodes and elements
  !
  nnodes = 0
  iflag = 0
  ierror = 0
  do while (iflag.eq.0) 
    read (8,'(a)', iostat=ierror1) tmp1
    ierror = ierror + ierror1
    if (index(tmp1, 'End Coordinates') > 0) then
      iflag = 1
    else
      nnodes = nnodes + 1
    endif
  enddo

  do i = 1,100
    read(8,*, iostat=ierror1) tmp1
    ierror = ierror + ierror1
    if (index(tmp1, 'Elements') > 0) exit
  end do

  iflag = 0
  nelems = 0

  do while (iflag.eq.0)
    read (8,'(a)', iostat=ierror1) tmp1
    ierror = ierror + ierror1
    if (index(tmp1, 'End Elements') > 0) then
      iflag = 1
    else
      nelems = nelems + 1 
    endif
  enddo
  rewind(8)
  if(ierror.gt.0) then
    ier = 2
    return
  endif

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
  

  read(8,'(a)', iostat=ierror1) tmp1
  read(8,'(a)', iostat=ierror2) tmp1

  if(ierror1+ierror2.gt.0) then
    ier = 2
    return
  endif
  
  ierror = 0
  do i=1,nnodes
    read(8,*, iostat=ierror1) node, Geometry1%Points(1,i), &
      Geometry1%Points(2,i), Geometry1%Points(3,i)
    ierror = ierror + ierror1
  enddo

  
  do i = 1,100
    read(8,*, iostat=ierror1) tmp1
    ierror = ierror + ierror1
    if (index(tmp1, 'Elements') > 0) exit
  end do

  do i=1,nelems
    read (8, *, iostat=ierror1) ielem, Geometry1%Tri(1,i), &
    Geometry1%Tri(2,i), Geometry1%Tri(3,i), Geometry1%Tri(4,i), &
    Geometry1%Tri(5,i), Geometry1%Tri(6,i)
    ierror = ierror + ierror1
  enddo

  if(ierror.gt.0) then
    ier = 2
    return
  endif


  close(8)
  

  return
end subroutine readgidmsh
!
!
!
!
!

subroutine readtri(Geometry1,filename, norder_skel, norder_smooth, ier)
  use ModType_Smooth_Surface
  implicit none

  !! This subroutine open a msh file and load the information in a
  !! variable of type Geometry

  type (Geometry), intent(inout) :: Geometry1     !! where the geometry will be loaded
  character(len=*), intent(in) :: filename         !! name of the msh file
  integer :: norder_smooth, norder_skel

  integer umio,i,m,N,j,aux1,aux2,aux3,ipointer
  integer :: ierror, nsk,nsf, ierror1, ierror2 
  integer :: ier

  ! set the flag for flat vs quadratic
  Geometry1%ifflat = 1

  
  ier = 0
  open(UNIT=8, FILE=trim(filename), STATUS='OLD', ACTION='READ', &
    IOSTAT=ierror)
  read(8,*, iostat=ierror1) m, N
  if(ierror+ierror1.gt.0) then
    ier = 2
    return
  endif

  nsk = (norder_skel+1)*(norder_skel+2)/2

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

  ierror = 0
  do j=1,m
    read(8,*, iostat=ierror1) Geometry1%Points(1,j), &
       Geometry1%Points(2,j), Geometry1%Points(3,j)
    ierror = ierror + ierror1
  enddo
  if(ierror.gt.0) then
    ier = 2
    return
  endif


  ipointer=m+1

  do j=1,N
    read(8,*, iostat=ierror1) aux1,aux2,aux3
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
    ierror = ierror + ierror1
  enddo
  if(ierror.gt.0) then
     ier = 2
     return
  endif
  close (8)

  return
end subroutine readtri




subroutine read_gmsh_v2(Geometry1, filename, norder_skel, norder_smooth, ier)
  use ModType_Smooth_Surface
  implicit none

  ! This subroutine opens a v2 gmsh file and load the information in a
  ! variable of type Geometry. Reads in all elements which are first
  ! or second order triangles, and first or second order quads
  !

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
  integer :: itype_tri3, itype_quad4
  integer :: itype_tri6, itype_quad8, itype_quad9
  integer :: inode1, inode2, inode3, inode4, npts, npts_use
  integer :: inode8(8) 
  integer :: nel_quad4, nel_quad9, nel_quad8, nel_tri3, nel_tri6
  integer :: ier

  ier = 0
  Geometry1%ifflat = 0

  iunit = 899

  open(UNIT=iunit, FILE=trim(filename), STATUS='OLD', ACTION='READ', IOSTAT=ierror)

  if(ierror.ne.0) then
    ier = 0
    return
  endif


  itype_tri3 = 2
  itype_tri6 = 9
  itype_quad4 = 4
  itype_quad9 = 10

!
!   Read the number of points
!
!

  do
    read(iunit, *, iostat=io) cline

    if (io .ne. 0) exit

    if (trim(cline) .eq. '$Nodes') then
      print *, 'Reading nodes . . . '
      read(iunit,*,iostat=io) numnodes
      if(io.ne.0) exit 
      print *, 'Number of nodes = ', numnodes
      print *
      
      allocate(xs(numnodes),ys(numnodes),zs(numnodes))
      do i = 1,numnodes
        read (iunit,*,iostat=io) ind, x, y, z
        xs(i) = x
        ys(i) = y
        zs(i) = z
      end do
      exit

    end if
  enddo

  ntri = 0
  npts = numnodes

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
        if(io.ne.0) exit 
        read (cline,*) ind, ielem, ntag

        if (ielem .eq. itype_tri3) then
          ntri = ntri + 1
          npts = npts + 3
        endif
        if (ielem .eq. itype_tri6) ntri = ntri + 1
        if (ielem .eq. itype_quad4) then
          ntri = ntri + 2
          npts = npts + 5
        endif
        if (ielem .eq. itype_quad8) then
          ntri = ntri + 2
          npts = npts + 1
        endif
        if (ielem .eq. itype_quad9) ntri = ntri + 2
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

  m = npts
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

  do j=1,numnodes
    Geometry1%Points(1,j) = xs(j)
    Geometry1%Points(2,j) = ys(j)
    Geometry1%Points(3,j) = zs(j)
  enddo

  npts_use = numnodes
  nel_tri3 = 3
  nel_tri6 = 6
  nel_quad4 = 4
  nel_quad9 = 9
  ntri = 0
  allocate(element(100))
  do
    read(iunit, *, iostat=io) cline
    if (io .ne. 0) exit

    if (trim(cline) .eq. '$Elements') then
      read(iunit,*) numelem
      do i = 1, numelem

        read(iunit, '(a)', iostat=io) cline
        read (cline,*) ind, ielem, ntag
        if (ielem.eq.itype_tri3) then
          nel = 3 + ntag + nel_tri3
          read(cline,*) element(1:nel)
          inode1 = element(3+ntag+1)
          inode2 = element(3+ntag+2)
          inode3 = element(3+ntag+3)

          ntri = ntri + 1
          Geometry1%Tri(1,ntri) = inode1
          Geometry1%Tri(2,ntri) = inode2
          Geometry1%Tri(3,ntri) = inode3
          Geometry1%Tri(4,ntri) = npts_use + 1
          Geometry1%Tri(5,ntri) = npts_use + 2
          Geometry1%Tri(6,ntri) = npts_use + 3

          npts_use = npts_use+1
          Geometry1%Points(1,npts_use) = (xs(inode1) + xs(inode2))/2
          Geometry1%Points(2,npts_use) = (ys(inode1) + ys(inode2))/2
          Geometry1%Points(3,npts_use) = (zs(inode1) + zs(inode2))/2
          
          npts_use = npts_use+1
          Geometry1%Points(1,npts_use) = (xs(inode2) + xs(inode3))/2
          Geometry1%Points(2,npts_use) = (ys(inode2) + ys(inode3))/2
          Geometry1%Points(3,npts_use) = (zs(inode2) + zs(inode3))/2
          
          npts_use = npts_use+1
          Geometry1%Points(1,npts_use) = (xs(inode3) + xs(inode1))/2
          Geometry1%Points(2,npts_use) = (ys(inode3) + ys(inode1))/2
          Geometry1%Points(3,npts_use) = (zs(inode3) + zs(inode1))/2
        endif

        if (ielem.eq.itype_quad4) then
          nel = 3 + ntag + nel_quad4
          read(cline,*) element(1:nel)
          inode1 = element(3+ntag+1)
          inode2 = element(3+ntag+2)
          inode3 = element(3+ntag+3)
          inode4 = element(3+ntag+4)


          ntri = ntri + 1
          Geometry1%Tri(1,ntri) = inode1
          Geometry1%Tri(2,ntri) = inode2
          Geometry1%Tri(3,ntri) = inode4
          Geometry1%Tri(4,ntri) = npts_use + 1
          Geometry1%Tri(5,ntri) = npts_use + 5
          Geometry1%Tri(6,ntri) = npts_use + 4

          ntri = ntri + 1
          Geometry1%Tri(1,ntri) = inode2
          Geometry1%Tri(2,ntri) = inode4
          Geometry1%Tri(3,ntri) = inode3
          Geometry1%Tri(4,ntri) = npts_use + 5
          Geometry1%Tri(5,ntri) = npts_use + 3
          Geometry1%Tri(6,ntri) = npts_use + 2


          npts_use = npts_use+1
          Geometry1%Points(1,npts_use) = (xs(inode1) + xs(inode2))/2
          Geometry1%Points(2,npts_use) = (ys(inode1) + ys(inode2))/2
          Geometry1%Points(3,npts_use) = (zs(inode1) + zs(inode2))/2
          
          npts_use = npts_use+1
          Geometry1%Points(1,npts_use) = (xs(inode2) + xs(inode3))/2
          Geometry1%Points(2,npts_use) = (ys(inode2) + ys(inode3))/2
          Geometry1%Points(3,npts_use) = (zs(inode2) + zs(inode3))/2
          
          npts_use = npts_use+1
          Geometry1%Points(1,npts_use) = (xs(inode3) + xs(inode4))/2
          Geometry1%Points(2,npts_use) = (ys(inode3) + ys(inode4))/2
          Geometry1%Points(3,npts_use) = (zs(inode3) + zs(inode4))/2

          npts_use = npts_use+1
          Geometry1%Points(1,npts_use) = (xs(inode4) + xs(inode1))/2
          Geometry1%Points(2,npts_use) = (ys(inode4) + ys(inode1))/2
          Geometry1%Points(3,npts_use) = (zs(inode4) + zs(inode1))/2

          npts_use = npts_use+1
          Geometry1%Points(1,npts_use) = (xs(inode1) + xs(inode2) + &
                                          xs(inode3) + xs(inode4))/4
          Geometry1%Points(2,npts_use) = (ys(inode1) + ys(inode2) + &
                                          ys(inode3) + ys(inode4))/4
          Geometry1%Points(3,npts_use) = (zs(inode1) + zs(inode2) + &
                                          zs(inode3) + zs(inode4))/4
          
        endif

        if (ielem.eq.itype_tri6) then
          nel = 3 + ntag + nel_tri6
          read(cline,*) element(1:nel)
          ntri = ntri + 1
          do j = 1,6
            Geometry1%Tri(j,ntri) = element(j+3+ntag)
          end do
        end if

        if (ielem.eq.itype_quad8) then
          nel = 3 + ntag + nel_quad8
          read(cline,*) element(1:nel)
          inode8(1:8) = element((3+ntag+1):(3+ntag+8))
          
          ntri = ntri + 1
          Geometry1%Tri(1,ntri) = element(3+ntag+1)
          Geometry1%Tri(2,ntri) = element(3+ntag+2)
          Geometry1%Tri(3,ntri) = element(3+ntag+4)
          Geometry1%Tri(4,ntri) = element(3+ntag+5)
          Geometry1%Tri(5,ntri) = npts_use + 1 
          Geometry1%Tri(6,ntri) = element(3+ntag+8)


          ntri = ntri + 1
          Geometry1%Tri(1,ntri) = element(3+ntag+3)
          Geometry1%Tri(2,ntri) = element(3+ntag+4)
          Geometry1%Tri(3,ntri) = element(3+ntag+2)
          Geometry1%Tri(4,ntri) = element(3+ntag+7)
          Geometry1%Tri(5,ntri) = npts_use + 1 
          Geometry1%Tri(6,ntri) = element(3+ntag+6)

          npts_use = npts_use + 1
          Geometry1%Points(1:3,npts_use) = 0
          do j=1,8
            Geometry1%Points(1,npts_use) = Geometry1%Points(1,npts_use) + &
              xs(inode8(j))
            Geometry1%Points(2,npts_use) = Geometry1%Points(1,npts_use) + &
              ys(inode8(j))
            Geometry1%Points(3,npts_use) = Geometry1%Points(1,npts_use) + &
              zs(inode8(j))
          enddo
          Geometry1%Points(1:3,npts_use) = Geometry1%Points(1:3,npts_use)/8

        endif
        

        if (ielem.eq.itype_quad9) then
          nel = 3 + ntag + nel_quad9
          read(cline,*) element(1:nel)
          ntri = ntri + 1
          Geometry1%Tri(1,ntri) = element(3+ntag+1)
          Geometry1%Tri(2,ntri) = element(3+ntag+2)
          Geometry1%Tri(3,ntri) = element(3+ntag+4)
          Geometry1%Tri(4,ntri) = element(3+ntag+5)
          Geometry1%Tri(5,ntri) = element(3+ntag+9)
          Geometry1%Tri(6,ntri) = element(3+ntag+8)


          ntri = ntri + 1
          Geometry1%Tri(1,ntri) = element(3+ntag+3)
          Geometry1%Tri(2,ntri) = element(3+ntag+4)
          Geometry1%Tri(3,ntri) = element(3+ntag+2)
          Geometry1%Tri(4,ntri) = element(3+ntag+7)
          Geometry1%Tri(5,ntri) = element(3+ntag+9)
          Geometry1%Tri(6,ntri) = element(3+ntag+6)

        endif
        
      end do

    end if
    
  end do
  deallocate(element)

  close(iunit)

  return
end subroutine read_gmsh_v2
!
!
!
!
!

subroutine read_gmsh_v4(Geometry1, filename, norder_skel, norder_smooth, ier)
  use ModType_Smooth_Surface
  implicit none

  ! This subroutine opens a v4 gmsh file and load the information in a
  ! variable of type Geometry. Reads in all elements which are first
  ! or second order triangles, and first or second order quads
  !

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
  integer :: itype_tri3, itype_quad4
  integer :: itype_tri6, itype_quad8, itype_quad9
  integer :: inode1, inode2, inode3, inode4, npts, npts_use
  integer :: inode8(8) 
  integer :: nel_quad4, nel_quad9, nel_quad8, nel_tri3, nel_tri6
  integer :: num_entity_blocks, mintag, maxtag
  integer :: ient, ienttag, iparam, iind, l, nelem, ielemtype
  integer, allocatable :: iindvec(:)

  integer ier


  Geometry1%ifflat = 0

  iunit = 899

  open(UNIT=iunit, FILE=trim(filename), STATUS='OLD', ACTION='READ', IOSTAT=ierror)
  if(ierror.ne.0) then
    ier = 2
    return
  endif


  itype_tri3 = 2
  itype_tri6 = 9
  itype_quad4 = 4
  itype_quad9 = 10

!
!   Read the number of points
!
!

  do
    read(iunit, *, iostat=io) cline

    if (io .ne. 0) exit

    if (trim(cline) .eq. '$Nodes') then
      print *, 'Reading nodes . . . '
      read(iunit,*) num_entity_blocks, numnodes, mintag, maxtag 
      print *, 'Number of nodes = ', numnodes
      print *

      allocate(xs(numnodes),ys(numnodes),zs(numnodes))
      allocate(iindvec(numnodes))
      do i=1, num_entity_blocks
         read(iunit, *, iostat = io) ient, ienttag, iparam, nnodes
         do j=1,nnodes
           read(iunit,*) iindvec(j)
         enddo
         do j=1,nnodes
           iind = iindvec(j)
           read(iunit,*) xs(iind), ys(iind), zs(iind)
         enddo
      enddo
      exit

    end if
  enddo

  ntri = 0
  npts = numnodes

  do
    read(iunit, *, iostat=io) cline
    if (io .ne. 0) exit

    if (trim(cline) .eq. '$Elements') then
      print *, 'Reading elements . . . '
      read(iunit,*) num_entity_blocks, numelem, mintag, maxtag 
      print *, 'Number of elements = ', numelem
      print *
      do i=1, num_entity_blocks
         read(iunit,*,iostat=io) ient, ienttag, ielemtype, nelem
         if(ielemtype.eq.itype_tri3) then
            ntri = ntri + nelem
            npts = npts + 3*nelem
         endif

         if(ielemtype.eq.itype_tri6) then
            ntri = ntri + nelem
         endif

         if(ielemtype.eq.itype_quad4) then
            ntri = ntri + 2*nelem
            npts = npts + 5*nelem
         endif

         if(ielemtype.eq.itype_quad8) then
            ntri = ntri + 2*nelem
            npts = npts + nelem
         endif

         if(ielemtype.eq.itype_quad9) then
            ntri = ntri + 2*nelem
         endif
         do j=1,nelem
           read(iunit,'(a)',iostat=io) cline
         enddo
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

  m = npts
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

  do j=1,numnodes
    Geometry1%Points(1,j) = xs(j)
    Geometry1%Points(2,j) = ys(j)
    Geometry1%Points(3,j) = zs(j)
  enddo

  npts_use = numnodes
  nel_tri3 = 3
  nel_tri6 = 6
  nel_quad4 = 4
  nel_quad9 = 9
  ntri = 0
  allocate(element(12))
  do
    read(iunit, *, iostat=io) cline
    if (io .ne. 0) exit

    if (trim(cline) .eq. '$Elements') then
      read(iunit,*) num_entity_blocks, numelem, mintag, maxtag 
      print *
      do i=1, num_entity_blocks
         read(iunit,*,iostat=io) ient, ienttag, ielemtype, nelem
         if(ielemtype.eq.itype_tri3) then
            nel = 1 + nel_tri3
            do j=1,nelem
              read(iunit,*) element(1:nel)
              inode1 = element(2)
              inode2 = element(3)
              inode3 = element(4)

              ntri = ntri + 1
              Geometry1%Tri(1,ntri) = inode1
              Geometry1%Tri(2,ntri) = inode2
              Geometry1%Tri(3,ntri) = inode3
              Geometry1%Tri(4,ntri) = npts_use + 1
              Geometry1%Tri(5,ntri) = npts_use + 2
              Geometry1%Tri(6,ntri) = npts_use + 3

              npts_use = npts_use+1
              Geometry1%Points(1,npts_use) = (xs(inode1) + xs(inode2))/2
              Geometry1%Points(2,npts_use) = (ys(inode1) + ys(inode2))/2
              Geometry1%Points(3,npts_use) = (zs(inode1) + zs(inode2))/2
          
              npts_use = npts_use+1
              Geometry1%Points(1,npts_use) = (xs(inode2) + xs(inode3))/2
              Geometry1%Points(2,npts_use) = (ys(inode2) + ys(inode3))/2
              Geometry1%Points(3,npts_use) = (zs(inode2) + zs(inode3))/2
          
              npts_use = npts_use+1
              Geometry1%Points(1,npts_use) = (xs(inode3) + xs(inode1))/2
              Geometry1%Points(2,npts_use) = (ys(inode3) + ys(inode1))/2
              Geometry1%Points(3,npts_use) = (zs(inode3) + zs(inode1))/2
            enddo
         endif

         if(ielemtype.eq.itype_tri6) then
           nel = 1 + nel_tri6
           do j=1,nelem
             read(iunit,*) element(1:nel)
             ntri = ntri + 1
             do l = 1,6
               Geometry1%Tri(l,ntri) = element(l+1)
             enddo
           enddo
         endif

         if(ielemtype.eq.itype_quad4) then
            nel = 1 + nel_quad4
            do j=1,nelem
              read(iunit,*) element(1:nel)
              inode1 = element(2)
              inode2 = element(3)
              inode3 = element(4)
              inode4 = element(5)

              ntri = ntri + 1
              Geometry1%Tri(1,ntri) = inode1
              Geometry1%Tri(2,ntri) = inode2
              Geometry1%Tri(3,ntri) = inode4
              Geometry1%Tri(4,ntri) = npts_use + 1
              Geometry1%Tri(5,ntri) = npts_use + 5
              Geometry1%Tri(6,ntri) = npts_use + 4

              ntri = ntri + 1
              Geometry1%Tri(1,ntri) = inode2
              Geometry1%Tri(2,ntri) = inode4
              Geometry1%Tri(3,ntri) = inode3
              Geometry1%Tri(4,ntri) = npts_use + 5
              Geometry1%Tri(5,ntri) = npts_use + 3
              Geometry1%Tri(6,ntri) = npts_use + 2


              npts_use = npts_use+1
              Geometry1%Points(1,npts_use) = (xs(inode1) + xs(inode2))/2
              Geometry1%Points(2,npts_use) = (ys(inode1) + ys(inode2))/2
              Geometry1%Points(3,npts_use) = (zs(inode1) + zs(inode2))/2
          
              npts_use = npts_use+1
              Geometry1%Points(1,npts_use) = (xs(inode2) + xs(inode3))/2
              Geometry1%Points(2,npts_use) = (ys(inode2) + ys(inode3))/2
              Geometry1%Points(3,npts_use) = (zs(inode2) + zs(inode3))/2
          
              npts_use = npts_use+1
              Geometry1%Points(1,npts_use) = (xs(inode3) + xs(inode4))/2
              Geometry1%Points(2,npts_use) = (ys(inode3) + ys(inode4))/2
              Geometry1%Points(3,npts_use) = (zs(inode3) + zs(inode4))/2

              npts_use = npts_use+1
              Geometry1%Points(1,npts_use) = (xs(inode4) + xs(inode1))/2
              Geometry1%Points(2,npts_use) = (ys(inode4) + ys(inode1))/2
              Geometry1%Points(3,npts_use) = (zs(inode4) + zs(inode1))/2

              npts_use = npts_use+1
              Geometry1%Points(1,npts_use) = (xs(inode1) + xs(inode2) + &
                                          xs(inode3) + xs(inode4))/4
              Geometry1%Points(2,npts_use) = (ys(inode1) + ys(inode2) + &
                                              ys(inode3) + ys(inode4))/4
              Geometry1%Points(3,npts_use) = (zs(inode1) + zs(inode2) + &
                                              zs(inode3) + zs(inode4))/4
            enddo
         endif

         if(ielemtype.eq.itype_quad8) then
            nel = 1 + nel_quad8
            do j=1,nelem
              read(iunit,*) element(1:nel)
              inode8(1:8) = element(2:9)
              
              ntri = ntri + 1
              Geometry1%Tri(1,ntri) = element(2)
              Geometry1%Tri(2,ntri) = element(3)
              Geometry1%Tri(3,ntri) = element(5)
              Geometry1%Tri(4,ntri) = element(6)
              Geometry1%Tri(5,ntri) = npts_use + 1 
              Geometry1%Tri(6,ntri) = element(9)


              ntri = ntri + 1
              Geometry1%Tri(1,ntri) = element(4)
              Geometry1%Tri(2,ntri) = element(5)
              Geometry1%Tri(3,ntri) = element(3)
              Geometry1%Tri(4,ntri) = element(8)
              Geometry1%Tri(5,ntri) = npts_use + 1 
              Geometry1%Tri(6,ntri) = element(7)

              npts_use = npts_use + 1
              Geometry1%Points(1:3,npts_use) = 0
              do l=1,8
                Geometry1%Points(1,npts_use) = Geometry1%Points(1,npts_use) + &
                  xs(inode8(l))
                Geometry1%Points(2,npts_use) = Geometry1%Points(1,npts_use) + &
                  ys(inode8(l))
                Geometry1%Points(3,npts_use) = Geometry1%Points(1,npts_use) + &
                  zs(inode8(l))
              enddo
              Geometry1%Points(1:3,npts_use) = Geometry1%Points(1:3,npts_use)/8
            enddo
         endif

         if(ielemtype.eq.itype_quad9) then
            nel = 1 + nel_quad9
            do j=1,nelem
              read(iunit,*) element(1:nel)
              ntri = ntri + 1
              Geometry1%Tri(1,ntri) = element(2)
              Geometry1%Tri(2,ntri) = element(3)
              Geometry1%Tri(3,ntri) = element(5)
              Geometry1%Tri(4,ntri) = element(6)
              Geometry1%Tri(5,ntri) = element(10)
              Geometry1%Tri(6,ntri) = element(9)


              ntri = ntri + 1
              Geometry1%Tri(1,ntri) = element(4)
              Geometry1%Tri(2,ntri) = element(5)
              Geometry1%Tri(3,ntri) = element(3)
              Geometry1%Tri(4,ntri) = element(8)
              Geometry1%Tri(5,ntri) = element(10)
              Geometry1%Tri(6,ntri) = element(7)
            enddo
         endif
      enddo

    end if
  end do
  deallocate(element)

  close(iunit)

  return
end subroutine read_gmsh_v4
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
!        * ifiletype = 4, for .msh gmsh v2
!        * ifiletype = 5, for .msh gmsh v4
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
  open(unit=33, file=trim(filename), status='old', iostat=io)
  if(io.ne.0) then
    ier = 1
    return

  endif

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
