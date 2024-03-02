subroutine plotskeletonvtk(Geometry1, filename)
  use Mod_Smooth_Surface
  implicit none

  type (Geometry) :: Geometry1
  character (len=*) filename

  integer :: umio,count1,count2,flag,n_order_sf, norder_smooth
  integer :: ierror, id, norder, nover, nsub, k, ntri, i, j, ictr
  integer :: ntot, ltot, npols7, info, iii, n, l, nnn, iw
  integer :: ifflat, is(10)
!  real (kind = 8) :: us(100000), vs(100000), ws(100000), dcond
!  real (kind = 8) :: uv1(10), uv2(10), uv3(10), uv(10), pols(100000)
!  real (kind = 8) :: xcoefs(10000), xrhs(10000)
!  real (kind = 8) :: ycoefs(10000), yrhs(10000)
!  real (kind = 8) :: zcoefs(10000), zrhs(10000)
  real (kind = 8) :: xval, yval, zval, pinv(1000000)

  real (kind = 8), allocatable :: xyzs(:,:,:), uvs(:,:,:)
  real (kind = 8), allocatable :: pmat(:,:), triout(:,:,:)

!  double precision :: umatr(100000), vmatr(100000)
  integer :: itype, npols
  integer ipatch,iunit1
  real *8 rr
  
  !
  ! This routien dumps out the skeleton patches (flat or quadratic)
  ! into vtk file,
  ! oversampling the triangles as necessary to show the smoothness
  !
  ! Input:
  !   Geometry1 - the structure containing all info
  !   filename - VTK ASCII filename, should end in .vtk
  !
  ! Output:
  !   the file 'filename' is created and contains vtk info
  !

  ifflat = Geometry1%ifflat

  call prinf('plotting skeleton, ifflat = *', ifflat, 1)

  if (ifflat .eq. 0) k = 6
  if (ifflat .eq. 1) k = 3

  !
  ! now dump out all the info needed for the triangles
  !
  ntri = Geometry1%ntri
  call prinf('in vtk plotter, original ntri = *', ntri, 1)
  
  allocate(triout(3,k,ntri))
  
  ! get vertex indices and construct triangles


  iunit1 = 877
  open(unit=iunit1,file=trim(filename),status='replace')
  write(iunit1,'(a)') "# vtk DataFile Version 3.0"
  write(iunit1,'(a)') "Skeleton"
  write(iunit1,'(a)') "ASCII"
  write(iunit1,'(a)') "DATASET UNSTRUCTURED_GRID"
  write(iunit1,'(a,i9,a)') "POINTS ", Geometry1%npoints, " float"
  do i=1,Geometry1%npoints
    write(iunit1,"(E11.5,2(2x,e11.5))") Geometry1%Points(1,i), &
      Geometry1%Points(2,i), Geometry1%Points(3,i)
  enddo

  write(iunit1,'(a,i9,i9)') "CELLS ", ntri, ntri*(k+1)

  if(ifflat.eq.0) then
    do i=1,ntri
      write(iunit1,'(a,i9,i9,i9,i9,i9,i9)') "6 ", Geometry1%Tri(1,i)-1, &
       Geometry1%Tri(2,i)-1,Geometry1%Tri(3,i)-1, &
       Geometry1%Tri(4,i)-1,Geometry1%Tri(5,i)-1, Geometry1%Tri(6,i)-1
    enddo
    write(iunit1,'(a,i9)') "CELL_TYPES ", ntri
    do ipatch = 1,ntri
      write(iunit1,'(a)') "22"
    end do

  elseif(ifflat.eq.1) then
    do i=1,ntri
      write(iunit1,'(a,i9,i9,i9)') "3 ", Geometry1%Tri(1,i)-1, &
       Geometry1%Tri(2,i)-1,Geometry1%Tri(3,i)-1
    enddo
    write(iunit1,'(a,i9)') "CELL_TYPES ", ntri
    do ipatch = 1,ntri
      write(iunit1,'(a)') "5"
    end do
  endif


  write(iunit1,'(a)') ""
  write(iunit1,'(a,i9)') "POINT_DATA ", Geometry1%npoints
  write(iunit1,'(a,i4)') "SCALARS normals float ", 3
  write(iunit1,'(a)') "LOOKUP_TABLE default"
  do i = 1,Geometry1%npoints
    write(iunit1,'(E11.5,2x,E11.5,2x,e11.5)') &
      Geometry1%Normal_Vert(1,i),&
      Geometry1%Normal_Vert(2,i),&
      Geometry1%Normal_Vert(3,i)
    rr = Geometry1%Normal_Vert(1,i)**2 + &
         Geometry1%Normal_Vert(2,i)**2 + &
         Geometry1%Normal_Vert(3,i)**2
  end do

  close(iunit1)



  return
end subroutine plotskeletonvtk





subroutine plotsmoothgeometryvtk(Geometry1, filename)
  use Mod_Smooth_Surface
  implicit none

  type (Geometry) :: Geometry1
  character (len=*) filename

  integer :: umio,count1,count2,flag,n_order_sf, norder_smooth
  integer :: ierror, id, norder, nover, nsub, k, ntri, i, j, ictr
  integer :: ntot, ltot, npols7, info, iii, n, l, nnn, iw
  real (kind = 8) :: dcond
  real (kind = 8) :: uv1(10), uv2(10), uv3(10), uv(10) 
  real (kind = 8), allocatable :: us(:),vs(:),ws(:),umatr(:,:),vmatr(:,:)
  real (kind = 8), allocatable :: uvtmp(:,:)
  real (kind = 8), allocatable :: xcoefs(:), xrhs(:)
  real (kind = 8), allocatable :: ycoefs(:), yrhs(:),pols(:)
  real (kind = 8), allocatable :: zcoefs(:), zrhs(:),pinv(:,:)
  real (kind = 8) :: xval, yval, zval

  real (kind = 8), allocatable :: xyzs(:,:,:), uvs(:,:,:)
  real (kind = 8), allocatable :: pmat(:,:), triout(:,:,:)

  integer :: itype, npols
  
  !
  ! This routien dumps out smoothed geometry into a vtk file,
  ! oversampling the triangles as necessary to show the smoothness
  !
  ! Input:
  !   Geometry1 - the structure containing all info
  !   filename - VTK ASCII filename, should end in .vtk
  !
  ! Output:
  !   the file 'filename' is created and contains vtk info
  !

  !id = 888
  !open(id, FILE=trim(filename),STATUS='REPLACE')

  norder_smooth = Geometry1%norder_smooth
  n_order_sf = Geometry1%n_order_sf

  !
  ! get the nodes here
  !

  norder = norder_smooth
  npols = (norder_smooth+1)*(norder_smooth+2)/2
  k = npols
  
  allocate(us(k),vs(k),umatr(k,k),vmatr(k,k),ws(k), uvtmp(2,k))
  allocate(xrhs(k),yrhs(k),zrhs(k),pinv(k,k))
  allocate(xcoefs(k),ycoefs(k),zcoefs(k),pols(k))
  call vioreanu_simplex_quad(norder_smooth, npols, uvtmp, &
      umatr, vmatr, ws)
  us = uvtmp(1,:)
  vs = uvtmp(2,:)


  if (n_order_sf .gt. 4**0) nover = 1
  if (n_order_sf .gt. 4**1) nover = 2
  if (n_order_sf .gt. 4**2) nover = 3
  if (n_order_sf .gt. 4**3) nover = 4
  if (n_order_sf .gt. 4**4) nover = 5

  nover = nover
  nsub = 4**nover

  !
  ! now dump out all the info needed for the triangles, compute xtri
  ! coefficients, and resample and plot
  !
  ntri = Geometry1%ntri
  call prinf('in vtk plotter, original ntri = *', ntri, 1)

  
  allocate(xyzs(3,k,ntri))


  ictr = 0
  do i = 1,ntri
    do j = 1,k
      ictr = ictr + 1
      xyzs(1,j,i) = Geometry1%S_smooth(1,ictr)
      xyzs(2,j,i) = Geometry1%S_smooth(2,ictr)
      xyzs(3,j,i) = Geometry1%S_smooth(3,ictr)
    end do
  end do

  
  allocate(uvs(2,3,nsub))
  uvs(1,1,1) = 0
  uvs(2,1,1) = 0
  uvs(1,2,1) = 1
  uvs(2,2,1) = 0
  uvs(1,3,1) = 0
  uvs(2,3,1) = 1


  !
  ! if necessary, recursively subdivide the triangle - first construct
  ! all the uv points
  !
  if (nover .gt. 0) then

    ntot = 1
    do i = 1,nover

      ltot = ntot

      do j = 1,ltot
        uv1(1) = uvs(1,1,j)
        uv1(2) = uvs(2,1,j)
        uv2(1) = uvs(1,2,j)
        uv2(2) = uvs(2,2,j)
        uv3(1) = uvs(1,3,j)
        uv3(2) = uvs(2,3,j)

        uvs(1,1,j) = uv1(1)
        uvs(2,1,j) = uv1(2)
        uvs(1,2,j) = (uv1(1) + uv2(1))/2
        uvs(2,2,j) = (uv1(2) + uv2(2))/2
        uvs(1,3,j) = (uv1(1) + uv3(1))/2
        uvs(2,3,j) = (uv1(2) + uv3(2))/2

        ntot = ntot + 1
        uvs(1,1,ntot) = (uv1(1) + uv2(1))/2
        uvs(2,1,ntot) = (uv1(2) + uv2(2))/2
        uvs(1,2,ntot) = uv2(1)
        uvs(2,2,ntot) = uv2(2)
        uvs(1,3,ntot) = (uv2(1) + uv3(1))/2
        uvs(2,3,ntot) = (uv2(2) + uv3(2))/2

        ntot = ntot + 1
        uvs(1,1,ntot) = (uv1(1) + uv3(1))/2
        uvs(2,1,ntot) = (uv1(2) + uv3(2))/2
        uvs(1,2,ntot) = (uv2(1) + uv3(1))/2
        uvs(2,2,ntot) = (uv2(2) + uv3(2))/2
        uvs(1,3,ntot) = uv3(1)
        uvs(2,3,ntot) = uv3(2)

        ntot = ntot + 1
        uvs(1,1,ntot) = (uv2(1) + uv3(1))/2
        uvs(2,1,ntot) = (uv2(2) + uv3(2))/2
        uvs(1,2,ntot) = (uv1(1) + uv3(1))/2
        uvs(2,2,ntot) = (uv1(2) + uv3(2))/2
        uvs(1,3,ntot) = (uv1(1) + uv2(1))/2
        uvs(2,3,ntot) = (uv1(2) + uv2(2))/2

      end do
    end do

  end if

  !call prinf('total triangles, ntot = *', ntot, 1)
  !call prin2('uvs = *', uvs, 6*ntot)

  !
  ! now compute the koornwinder expansion of the triangle
  !
  allocate(pmat(k,k))

  !call prin2('us = *', us, k)
  !call prin2('vs = *', vs, k)
  !print *
  !print *

  do i = 1,k
    uv(1) = us(i)
    uv(2) = vs(i)
    call koorn_pols(uv, norder, npols, pols)
    !call prin2('uv = *', uv, 2)
    !call prinf('npols = *', npols, 1)
    if (npols .ne. n_order_sf) then
      call prinf('npols = *', npols, 1)
      call prinf('n_order_sf = *', n_order_sf, 1)
      stop
    end if
    
    !call prin2('pols = *', pols, k)
    !stop
    do j = 1,npols
      pmat(i,j) = pols(j)
    end do
  end do


  call dinverse(npols, pmat, info, pinv)
  !call prinf('after inverse, info = *', info, 1)
  
  !
  ! loop over each triangle, solve for each of the sets of
  ! coefficients, and then evaluate the subsampled triangles
  !

  allocate(triout(3,3,nsub*ntri))
  nnn = 0

  do i = 1,ntri

    do j = 1,k
      xrhs(j) = xyzs(1,j,i)
      yrhs(j) = xyzs(2,j,i)
      zrhs(j) = xyzs(3,j,i)
    end do

    call dmatvec(k, k, pinv, xrhs, xcoefs)
    call dmatvec(k, k, pinv, yrhs, ycoefs)
    call dmatvec(k, k, pinv, zrhs, zcoefs)

    !
    ! now evaluate the new triangle nodes
    !
    do j = 1,nsub
      
      nnn = nnn + 1
      
      do iii = 1,3
        uv(1) = uvs(1,iii,j)
        uv(2) = uvs(2,iii,j)
        call koorn_pols(uv, norder, npols, pols)
        xval = 0
        yval = 0
        zval = 0
        do l = 1,k
          xval = xval + xcoefs(l)*pols(l)
          yval = yval + ycoefs(l)*pols(l)
          zval = zval + zcoefs(l)*pols(l)
        end do
        
        triout(1,iii,nnn) = xval
        triout(2,iii,nnn) = yval
        triout(3,iii,nnn) = zval
        
      end do

      !call prin2('tri = *', triout(1,1,nnn), 9)
      
    end do
    
    
  end do

  call prinf('num triangles plotted = *', nnn, 1)
  
  
  call xtri_vtk_flat(nnn, triout, 'smoothed geometry', filename)

  !close (id)
  return
end subroutine plotsmoothgeometryvtk





subroutine xtri_vtk_flat(ntri, xtri1s, title, filename)
  implicit real *8 (a-h,o-z)
  real *8 :: xtri1s(3,3,ntri)
  character(len=*) :: title, filename

  character(len=1024) :: dataname, valsname, imgname
  character(len=1024) :: trisname, vecsname, centname
  character(len=12) :: fmt, fmt3, fmt4
  character(len=25) :: fmt2

  !
  ! This routine plots a sequence of FLAT TRIANGLES with surface
  ! color vals.
  !
  ! Input:
  !   ntri - number of flat triangles
  !   xtri1s - full triangle information
  !
  ! Output:
  !   files which can be executed in matlab to plot the surface
  !
  !

  !
  ! write the vtk plotting script
  !
  iunit1 = 877
  open(unit = iunit1, file=trim(filename), status='replace')

  write(iunit1,'(a)') "# vtk DataFile Version 3.0"
  write(iunit1,'(a)') "vtk output"
  write(iunit1,'(a)') "ASCII"
  !write(iunit1,'(a)') "DATASET POLYDATA"
  write(iunit1,'(a)') "DATASET UNSTRUCTURED_GRID"
  write(iunit1,'(a,i8,a)') "POINTS ", ntri*3, " float"

  fmt2 = "(E11.5,2X,E11.5,2X,E11.5)"
  do i = 1,ntri
    write(iunit1,fmt2) xtri1s(1,1,i), xtri1s(2,1,i), xtri1s(3,1,i)
    write(iunit1,fmt2) xtri1s(1,2,i), xtri1s(2,2,i), xtri1s(3,2,i)
    write(iunit1,fmt2) xtri1s(1,3,i), xtri1s(2,3,i), xtri1s(3,3,i)
  end do


  write(iunit1,'(a,i9,i9)') "CELLS ", ntri, ntri*4

  do i = 1,ntri
    i1 = 3*(i-1) + 1
    write(iunit1,'(a,i9,i9,i9)') "3 ", i1-1, i1, i1+1
  end do

  write(iunit1,'(a,i9)') "CELL_TYPES ", ntri
  do i = 1,ntri
    write(iunit1,'(a)') "5"
  end do

  write(iunit1,'(a,i8)') "POINT_DATA ", ntri*3
  write(iunit1,'(a)') "SCALARS scalars float 1"
  write(iunit1,'(a)') "LOOKUP_TABLE default"
  do i = 1,ntri
    do j = 1,3
      write(iunit1,'(E11.5)') xtri1s(3,j,i)
    end do
  end do



  write(iunit1,'(a)') ""
  write(iunit1,'(a)') ""
  write(iunit1,'(a,i9)') "CELL_DATA ", ntri
  write(iunit1,'(a)') "SCALARS scalars float 1"
  write(iunit1,'(a)') "LOOKUP_TABLE default"
  do i = 1,ntri
    write(iunit1,'(E13.5)') (xtri1s(3,1,i) + &
        xtri1s(3,2,i) + xtri1s(3,3,i))/3
  end do

  close(iunit1)

  return
end subroutine xtri_vtk_flat





subroutine xtri_vtk_quadratic(ntri, xtri1s, title, filename)
  implicit real *8 (a-h,o-z)
  real *8 :: xtri1s(3,6,ntri)
  character(len=*) :: title, filename

  character(len=1024) :: dataname, valsname, imgname
  character(len=1024) :: trisname, vecsname, centname
  character(len=12) :: fmt, fmt3, fmt4
  character(len=25) :: fmt2

  !
  ! This routine plots a sequence of FLAT TRIANGLES with surface
  ! color vals.
  !
  ! Input:
  !   ntri - number of flat triangles
  !   xtri1s - full triangle information
  !
  ! Output:
  !   files which can be executed in matlab to plot the surface
  !
  !

  !
  ! write the vtk plotting script
  !
  iunit1 = 877
  open(unit = iunit1, file=trim(filename), status='replace')

  write(iunit1,'(a)') "# vtk DataFile Version 3.0"
  write(iunit1,'(a)') "vtk output"
  write(iunit1,'(a)') "ASCII"
  !write(iunit1,'(a)') "DATASET POLYDATA"
  write(iunit1,'(a)') "DATASET UNSTRUCTURED_GRID"
  write(iunit1,'(a,i8,a)') "POINTS ", ntri*6, " float"

  fmt2 = "(E11.5,2X,E11.5,2X,E11.5)"
  do i = 1,ntri
    do j = 1,6
      write(iunit1,fmt2) xtri1s(1,j,i), xtri1s(2,j,i), xtri1s(3,j,i)
    end do
  end do


  write(iunit1,'(a,i8,i8)') "CELLS ", ntri, ntri*7

  do i = 1,ntri
    i1 = 6*(i-1) + 1
    write(iunit1,'(a,i8,i8,i8,i8,i8,i8)') "6 ", i1-1, i1, i1+1, &
        i1+2, i1+3, i1+4
  end do

  write(iunit1,'(a,i8)') "CELL_TYPES ", ntri
  do i = 1,ntri
    write(iunit1,'(a)') "22"
  end do

  write(iunit1,'(a,i8)') "POINT_DATA ", ntri*6
  write(iunit1,'(a)') "SCALARS scalars float 1"
  write(iunit1,'(a)') "LOOKUP_TABLE default"
  do i = 1,ntri
    do j = 1,6
      write(iunit1,'(E11.5)') xtri1s(3,j,i)
    end do
  end do

  write(iunit1,'(a)') ""
  write(iunit1,'(a)') ""
  write(iunit1,'(a,i8)') "CELL_DATA ", ntri
  write(iunit1,'(a)') "SCALARS scalars float 1"
  write(iunit1,'(a)') "LOOKUP_TABLE default"
  do i = 1,ntri
    write(iunit1,'(E13.5)') (xtri1s(3,1,i) + &
        xtri1s(3,2,i) + xtri1s(3,3,i))/3
  end do

  close(iunit1)

  return
end subroutine xtri_vtk_quadratic






subroutine plot1D(x,y,N,nombre)
implicit none

!List of calling arguments
character (len=100), intent(in) :: nombre
integer, intent(in) :: N
real ( kind = 8 ), intent(in) :: x(N),y(N)

!List of local variables
!character (len=100) nombre
integer umio,count,flag
integer :: ierror
!    nombre='fortran_plot'
    open(8, FILE=trim(nombre),STATUS='REPLACE')

!    open(UNIT=8, FILE=nombre, STATUS='OLD', ACTION='REPLACE', IOSTAT=ierror)
    flag=1
    write(8,*) flag
    write(8,*) N
    do count=1,N
        write(8,*) x(count)
    enddo
    do count=1,N
        write(8,*) y(count)
    enddo
    close (8)
return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine plot_curve_3D(x,y,z,N,nombre)
implicit none

!List of calling arguments
character (len=100), intent(in) :: nombre
integer, intent(in) :: N
real ( kind = 8 ), intent(in) :: x(N),y(N),z(N)

!List of local variables
!character (len=100) nombre
integer umio,count,flag
integer :: ierror
!    nombre='fortran_plot'
    open(8, FILE=trim(nombre),STATUS='REPLACE')

!    open(UNIT=8, FILE=nombre, STATUS='OLD', ACTION='REPLACE', IOSTAT=ierror)
    flag=3
    write(8,*) flag
    write(8,*) N
    do count=1,N
        write(8,*) x(count)
    enddo
    do count=1,N
        write(8,*) y(count)
    enddo
    do count=1,N
        write(8,*) z(count)
    enddo
    close (8)
return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plot2D(x,y,F,N,M,nombre)
implicit none

!List of calling arguments
character (len=100), intent(in) :: nombre
integer, intent(in) :: N,M
real ( kind = 8 ), intent(in) :: x(M,N),y(M,N),F(M,N)

!List of local variables
!character (len=100) nombre
integer umio,count1,count2,flag
integer :: ierror
!    nombre='fortran_plot'
    open(8, FILE=trim(nombre),STATUS='REPLACE')
    flag=2
    write(8,*) flag
    write(8,*) M
    write(8,*) N
    do count1=1,N
        do count2=1,M
            write(8,*) x(count2,count1)
        enddo
    enddo
    do count1=1,N
        do count2=1,M
            write(8,*) y(count2,count1)
        enddo
    enddo
    do count1=1,N
        do count2=1,M
            write(8,*) F(count2,count1)
        enddo
    enddo
    close (8)
return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plot2D_v2(F,N,M,nombre)
implicit none

!List of calling arguments
character (len=100), intent(in) :: nombre
integer, intent(in) :: N,M
real ( kind = 8 ), intent(in) :: F(M,N)

!List of local variables
!character (len=100) nombre
integer umio,count1,count2,flag
integer :: ierror
!    nombre='fortran_plot'
    open(8, FILE=trim(nombre),STATUS='REPLACE')
    flag=4
    write(8,*) flag
    write(8,*) M
    write(8,*) N
    do count1=1,N
        do count2=1,M
            write(8,*) F(count2,count1)
        enddo
    enddo
    close (8)
return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plot2D_v3(x,y,z,F,N,M,nombre)
implicit none

!List of calling arguments
character (len=100), intent(in) :: nombre
integer, intent(in) :: N,M
real ( kind = 8 ), intent(in) :: x(M,N),y(M,N),z(M,N),F(M,N)

!List of local variables
integer umio,count1,count2,flag
integer :: ierror
!    nombre='fortran_plot'
    open(8, FILE=trim(nombre),STATUS='REPLACE')
    flag=5
    write(8,*) flag
    write(8,*) M
    write(8,*) N
    do count1=1,N
        do count2=1,M
            write(8,*) x(count2,count1)
        enddo
    enddo
    do count1=1,N
        do count2=1,M
            write(8,*) y(count2,count1)
        enddo
    enddo
    do count1=1,N
        do count2=1,M
            write(8,*) z(count2,count1)
        enddo
    enddo
    do count1=1,N
        do count2=1,M
            write(8,*) F(count2,count1)
        enddo
    enddo
    close (8)
return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plot_tree(W_boxes,num_box,Pts,n_pts,nombre)
implicit none

!List of calling arguments
character (len=100), intent(in) :: nombre
integer, intent(in) :: num_box,n_pts
real ( kind = 8 ), intent(in) :: W_boxes(num_box*4),Pts(3,n_pts)

!List of local variables
integer umio,count1,count2,flag
integer :: ierror
!    nombre='fortran_plot'
    open(8, FILE=nombre,STATUS='REPLACE')
    flag=6
    write(8,*) flag
    write(8,*) num_box
    write(8,*) n_pts
    do count1=1,n_pts
        write(8,*) Pts(1,count1)
    enddo
    do count1=1,n_pts
        write(8,*) Pts(2,count1)
    enddo
    do count1=1,n_pts
        write(8,*) Pts(3,count1)
    enddo

    do count1=1,num_box
        write(8,*) W_boxes(4*(count1-1)+1)
        write(8,*) W_boxes(4*(count1-1)+2)
        write(8,*) W_boxes(4*(count1-1)+3)
        write(8,*) W_boxes(4*(count1-1)+4)
    enddo
    close (8)
return
end
