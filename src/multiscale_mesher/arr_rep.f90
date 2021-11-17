subroutine read_tri(norder_smooth, tol, filename, v_coord, tri_ind, tri_ind_se, vt_ind, vt_ind_se)
    implicit none
    double precision tol
    integer norder_smooth
    character(len=100), intent(in) :: filename         !! name of the msh file

    integer i,m,n,j,aux1,aux2,aux3,v_i
    integer, allocatable :: tri_ind(:,:), tri_ind_se(:,:)
    integer, allocatable :: vt_ind(:), vt_ind_se(:,:)
    integer, allocatable :: v_count(:), v_offset(:), v_icount(:)
    double precision, allocatable :: v_coord(:,:)

    open(UNIT=8, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ierror)
    read(8,*) m, n


    print *
    print *
    write (6,*) 'loading file ', trim(filename)
    write (13,*) 'loading file ', trim(filename)
    call prinf('npoints = *', m, 1)
    call prinf('ntri = *', n, 1)


    if (allocated(v_coord)) then
        deallocate(v_coord)
    endif
    if (allocated(tri_ind)) then
        deallocate(tri_ind)
    endif
    if (allocated(tri_ind_se)) then
        deallocate(tri_ind_se)
    endif
    if (allocated(vt_ind)) then
        deallocate(vt_ind)
    endif
    if (allocated(vt_ind_se)) then
        deallocate(vt_ind_se)
    endif

    allocate(v_coord(3,m))
    allocate(tri_ind(3,n))
    allocate(tri_ind_se(2,n))
    allocate(vt_ind_se(2,n))
    allocate(v_count(m))
    allocate(v_off(m))

    do j=1,m
        read(8,*) v_coord(1,j),v_coord(2,j),v_coord(3,j)
    enddo

    do j=1,n
        ! need to change tri_ind to one dimensional array
        read(8,*) tri_ind(1,j),tri_ind(2,j),tri_ind(3,j)
        tri_ind_se(1,j) = 3*(j-1) + 1
        tri_ind_se(2,j) = 3*j
    enddo
    close (8)

    do j=1,m
        v_count(j) = 0
        v_off(j) = 0
        v_icount(j) = 0
    enddo

    do j=1,n
        v_count(tri_ind(1,j)) = v_count(tri_ind(1,j)) + 1
        v_count(tri_ind(2,j)) = v_count(tri_ind(2,j)) + 1
        v_count(tri_ind(3,j)) = v_count(tri_ind(3,j)) + 1
    enddo
    
    do j=2,n
       v_off(j) = v_off(j-1) + v_count(j-1)
    enddo

    do j=1,m
        vt_ind_se(1,j) = v_off(j) + 1
        vt_ind_se(2,j) = vt_ind_se(1,j) + v_count(j)
    enddo

    allocate(vt_ind(vt_ind_se(m)))

    do j=1,n
    
        v_i = tri_ind(1,j)
        vt_ind(v_off(v_i) + v_icount(v_i) + 1) = tri_ind(1,j)
        v_icount(v_i) = v_icount(v_i) + 1

        v_i = tri_ind(2,j)
        vt_ind(v_off(v_i) + v_icount(v_i) + 1) = tri_ind(1,j)
        v_icount(v_i) = v_icount(v_i) + 1

        v_i = tri_ind(3,j)
        vt_ind(v_off(v_i) + v_icount(v_i) + 1) = tri_ind(1,j)
        v_icount(v_i) = v_icount(v_i) + 1

    enddo

end subroutine read_tri

subroutine compute_v_quantities(nv, nt, v_coord, tri_ind, tri_ind_se, vt_ind, vt_ind_se, tri_normal, v_normal, tri_theta)
    integer nv, nt
    double precision, allocatable :: v_coord(3,:)
    integer, allocatable :: tri_ind(3,:)
    integer, allocatable :: tri_ind_se(2,:)
    integer, allocatable :: vt_ind(:)
    integer, allocatable :: vt_ind_se(2,:)
    double precision, allocatable :: tri_normal(3,nt)
    double precision, allocatable :: v_normal(3,nv)
    double precision, allocatable :: tri_theta(3,nt)
    integer i,j,ai,bi,ci,ti
    double precision a,b,c,v_tmp(3),t1(3),t2(3),denom
    double precision theta,theta_sum

    ! theta and unit normal for each tri
    do i = 1,nt
        ai = tri_ind(1,i)
        bi = tri_ind(2,i)
        ci = tri_ind(3,i)
        
        a = 0.0d0
        do j = 1,3
            v_tmp(j) = v_coord(j,bi) - v_coord(j,ci)
            a = a + v_tmp(j)**2
        enddo
        a = sqrt(a)

        b = 0.0d0
        do j = 1,3
            v_tmp(j) = v_coord(j,ai) - v_coord(j,ci)
            b = b + v_tmp(j)**2
        enddo
        b = sqrt(b)

        c = 0.0d0
        do j = 1,3
            v_tmp(j) = v_coord(j,ai) - v_coord(j,bi)
            c = c + v_tmp(j)**2
        enddo
        c = sqrt(c)

        tri_theta(1,i) = acos((b**2+c**2-a**2)/(2*b*c))
        tri_theta(2,i) = acos((a**2+c**2-b**2)/(2*a*c))
        tri_theta(3,i) = acos((a**2+b**2-c**2)/(2*a*b))

        do j=1,3
            t1(j) = v_coord(j,bi) - v_coord(j,ai)
            t2(j) = v_coord(j,ci) - v_coord(j,ai)
        enddo
        v_tmp(1) = t1(2)*t2(3) - t1(3)*t2(2)
        v_tmp(2) = -t1(1)*t2(3) + t1(3)*t2(1)
        v_tmp(3) = t1(1)*t2(2) - t1(2)*t2(1)
        denom = sqrt(v_tmp(1)**2+v_tmp(2)**2+v_tmp(3)**2)
        do j=1,3
            tri_normal(j,i) = v_tmp(j)/denom
        enddo
        
    enddo

    do i = 1,nv
        theta_sum = 0.0d0
        do k = 1,3
            v_tmp(k) = 0.0d0
        enddo
        do j = vt_ind_se(1,i), vt_ind_se(2,i)
            ti = vt_ind(j)
            if(i==tri_ind(1,ti)) then
                theta = tri_theta(1,ti)
            endif
            if(i==tri_ind(2,ti)) then
                theta = tri_theta(2,ti)
            endif
            if(i==tri_ind(3,ti)) then
                theta = tri_theta(3,ti)
            endif
            theta_sum = theta_sum + theta
            do k = 1,3
                v_tmp(k) = theta*tri_normal(k,ti)
            enddo
        enddo
        do k = 1,3
            v_normal(k,i) = v_tmp(k)/theta_sum
        enddo
    enddo

end subroutine compute_v_quantities

subroutine read_q_gmsh(Geometry1, filename, norder_skel, norder_smooth)
  use ModType_Smooth_Surface
  implicit none

  !! This subroutine opens a v2 gmsh file and load the information in a
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
  character(len=100), intent(in) :: filename         !! name of the msh file
  character(len=100) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7
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


  Geometry1%ifflat = 0

  iunit = 899

  open(UNIT=iunit, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ierror)

  write (6,*) 'loading file ', trim(filename)
  write (13,*) 'loading file ', trim(filename)
  print *

  korder = 2
  kpols = (korder+1)*(korder+2)/2

  itype = -1
  if (korder .eq. 1) itype = 2
  if (korder .eq. 2) itype = 9
  if (korder .eq. 3) itype = 21
  if (korder .eq. 4) itype = 23
  if (korder .eq. 5) itype = 25
  if (korder .eq. 6) itype = -1






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

    end if

    if (trim(cline) .eq. '$Elements') then
      print *, 'Reading elements . . . '
      read(iunit,*) numelem
      print *, 'Number of elements = ', numelem
      print *

      nel = (korder+1)*(korder+2)/2
      allocate(elements(nel,numelem))

      ntri = 0
      do i = 1,numelem

        read(iunit, '(a)', iostat=io) cline
        read (cline,*) ind, ielem, ntag

        if (ielem .eq. itype) then

          ntri = ntri + 1

          lll= 1+1+1+ntag+nel
          allocate(element(lll))
          read(cline,*) element

          do j = 1,kpols
            elements(j,ntri) = element(j+3+ntag)
          end do
          deallocate(element)

        end if

      end do

    end if

  end do



  n = ntri
  call prinf('ntri = *', n, 1)

  nsk = (norder_skel+1)*(norder_skel+2)/2

  call prinf('num points on skeleton mesh = *', nsk*n, 1)

  Geometry1%norder_skel = norder_skel
  Geometry1%nskel = nsk


  Geometry1%norder_smooth = norder_smooth
  nsf = (norder_smooth+1)*(norder_smooth+2)/2
  Geometry1%nsmooth = nsf

  call prinf('num points on smooth mesh = *', nsf*n, 1)

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

  do j=1,N
    do i = 1,6
      Geometry1%Tri(i,j) = elements(i,j)
    end do
  enddo


  close(iunit)

  return
end subroutine read_q_gmsh
