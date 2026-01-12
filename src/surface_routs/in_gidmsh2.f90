!
!  This file contains the input routines for reading in 
!  Gmsh v2 file on input, and converting the input data structure
!  to a Rokhlin vioreanu discretization of an equivalent order.
!
!  The file currently supports reading low order triangular elements
!  with 3,6,10, or 15 nodes
!
!  The file has 2 user callable routines
!     open_gmsh2_geometry_mem - estimates the memory reuirements
!     open_gmsh2_geometry - initializes various arrays
!
!
!
!
subroutine open_gidmsh2_geometry(filename,npatches,norders,ixyzs, &
  iptype,npts,srcvals,srccoefs,wts)
!
! This subroutine computes the srcvals, srccoefs, norders, ixyzs,
! iptype, and wts array from a gmsh v2 file.
!
! The file reads in triangles discretized using 1,2...5 
! gmsh nodes, computes koornwinder expansion coefficients
! for the xyz coordinates from the coordinates provided by the gmsh
! file, resamples the koornwinder expansion at the same order
! Vioreanu Rokhlin discretization nodes, and obtains the derivative
! information on each patch via spectral diffrentiation.
!
!  Input arguments:
!    - filename: character(len=*)
!         file name of gmshv2 file
!    - npatches: integer *8
!         number of relevant triangles in gmsh file. Should be computed
!         using a call to open_gmshv2_geometry_mem
!    - npts: integer *8
!         total number of discretization points. Should be computed
!         using a call to open_gmshv2_geometry_mem
!
!  Output arguments:
!    
!    - norders: integer *8(npatches)
!        order of discretization of the patches
!    - ixyzs: integer *8(npatches+1)
!        starting location of points on patch i
!    - iptype: integer *8(npatches)
!        type of patch
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!    - wts: double precision (npts)
!        quadrature weights for integrating smooth functions
!        at discretization nodes
!

  implicit none
  character (len=*), intent(in) :: filename
  integer ( kind = 8 ), intent(in) :: npts,npatches
  integer ( kind = 8 ), intent(out) :: norders(npatches),ixyzs(npatches+1)
  integer ( kind = 8 ), intent(out) :: iptype(npatches)
  real ( kind = 8 ), intent(out) :: srcvals(12,npts),srccoefs(9,npts)
  real ( kind = 8 ), intent(out) :: wts(npts)

  ! List of local variables
  real ( kind = 8 ), allocatable :: umatall(:,:,:),uv0all(:,:,:)
  real ( kind = 8 ), allocatable :: uvvrall(:,:,:),wvrall(:,:)
  real ( kind = 8 ), allocatable :: umatvrall(:,:,:)
  real ( kind = 8 ), allocatable :: umat(:,:),vmat(:,:),pols(:)
  real ( kind = 8 ), allocatable :: pmatall(:,:,:),dmatall(:,:,:,:)

  real ( kind = 8 ), allocatable :: xyzs(:,:)
  integer ( kind = 8 ), allocatable :: element(:)
  integer *8 n,nomax,i,ind,info,inode,ipt,istart,itri,j,l,m
  integer *8 norder,npoints,iunit,nmax,npt,ntri,npols,ierror
  real *8 rr

 
!
!  Gather all node and vals2coefs (umats), and coefs2vals (pmat,dmat,vmat) 
!  matrices
!
 
  nomax = 2
  nmax = (nomax+1)*(nomax+2)/2
  allocate(uv0all(2,nmax,nomax))
  call load_uv_gmsh2_all(nomax,nmax,uv0all)
  allocate(pols(nmax))
  allocate(umatall(nmax,nmax,nomax),umatvrall(nmax,nmax,nomax))
  allocate(uvvrall(2,nmax,nomax),wvrall(nmax,nomax))
  call load_uvw_vr_all(nomax,nmax,uvvrall,wvrall)

  do i=1,nomax
    do j=1,nmax
      do l=1,nmax
        umatall(l,j,i) = 0
        umatvrall(l,j,i) = 0
      enddo
    enddo
  enddo

  do i=1,nomax
    norder = i
    npols = (norder+1)*(norder+2)/2
    allocate(umat(npols,npols),vmat(npols,npols))
    do j=1,npols
      call koorn_pols(uv0all(1,j,i),norder,npols,pols)
      do l=1,npols
        vmat(j,l) = pols(l)
      enddo
    enddo
    call dinverse(npols,vmat,info,umat)
    do j=1,npols
      do l=1,npols
        umatall(l,j,i) = umat(l,j)
      enddo
    enddo

    do j=1,npols
      call koorn_pols(uvvrall(1,j,i),norder,npols,pols)
      do l=1,npols
        vmat(j,l) = pols(l)
      enddo
    enddo
    call dinverse(npols,vmat,info,umat)
    do j=1,npols
      do l=1,npols
        umatvrall(l,j,i) = umat(l,j)
      enddo
    enddo
    deallocate(umat,vmat)
  enddo

  allocate(pmatall(nmax,nmax,nomax),dmatall(2,nmax,nmax,nomax))

  do i=1,nomax
    do j=1,nmax
      do l=1,nmax
        pmatall(l,j,i) = 0
        dmatall(1,l,j,i) = 0
        dmatall(2,l,j,i) = 0
      enddo
    enddo
    npols = (i+1)*(i+2)/2
    do j=1,npols
      call koorn_ders(uvvrall(1,j,i),i,npols,pmatall(1,j,i), &
        dmatall(1,1,j,i))
    enddo
  enddo

  iunit = 100
  !open(unit=iunit, file='../../xtri/geometry/cubetorus_nodes.txt', &
  !    status='old', action='read')

  open(unit=iunit, file=trim(filename), status='old', action='read')

  ixyzs(1) = 1
  read(iunit,*) ind,ind,ind,npt,ind
  allocate(xyzs(3,npt))
  do i=1,npt
    read(iunit,*) ind,ind,ind,ind,ind,ind,ind,ind
    read(iunit,*) xyzs(1,i),xyzs(2,i),xyzs(3,i)
  enddo
  read(iunit,*)
  norder = 2
  npols = (norder+1)*(norder+2)/2
  allocate(element(npols))
  do itri=1,npatches
    read(iunit,*) ind,ind,ind,ind,ind,ind,ind,ind
    read(iunit,*) element
    ixyzs(itri) = (itri-1)*npols + 1
!
!  Temporarily store xyz coordinates for the triangle in 
!  srcvals (even though it is still gmsh nodes) srcvals 
!  will be overriden later
!
    do j = 1,npols
      inode = element(j)
      srcvals(1:3,ixyzs(itri)+j-1) = xyzs(1:3,inode)
    end do
    istart = ixyzs(itri)
!
!  Form koornwinder expansion coeffs from xyzs info
!
    do l=1,npols
      srccoefs(1:3,istart+l-1) = 0
      do m=1,npols
        srccoefs(1:3,istart+l-1) = srccoefs(1:3,istart+l-1) + &
            umatall(l,m,norder)*srcvals(1:3,istart+m-1)
      enddo
    enddo
!
!
!  Evaluate koornwinder expansion coeffs at vioreanu rokhlin nodes
!

    do l=1,npols
      srcvals(1:9,istart+l-1) = 0
      do m=1,npols
        srcvals(1:3,istart+l-1) = srcvals(1:3,istart+l-1) + &
            pmatall(m,l,norder)*srccoefs(1:3,istart+m-1)
        srcvals(4:6,istart+l-1) = srcvals(4:6,istart+l-1) + &
            dmatall(1,m,l,norder)*srccoefs(1:3,istart+m-1)
        srcvals(7:9,istart+l-1) = srcvals(7:9,istart+l-1) + &
            dmatall(2,m,l,norder)*srccoefs(1:3,istart+m-1)
      enddo
    enddo
    
    do l=1,npols
      srccoefs(4:9,istart+l-1) = 0
      do m=1,npols
        srccoefs(4:9,istart+l-1) = srccoefs(4:9,istart+l-1) + &
           umatvrall(l,m,norder)*srcvals(4:9,istart+m-1)
      enddo
    enddo

    do l=1,npols
      ipt = istart + l-1
      call cross_prod3d(srcvals(4:6,ipt),srcvals(7:9,ipt), &
           srcvals(10:12,ipt))
      rr = sqrt(srcvals(10,ipt)**2 + srcvals(11,ipt)**2 + &
           srcvals(12,ipt)**2)
      wts(ipt) = rr*wvrall(l,norder)
      srcvals(10:12,ipt) = srcvals(10:12,ipt)/rr
    enddo
    norders(itri) = norder
    iptype(itri) = 1
  end do
  ixyzs(npatches+1) = npts+1

  close(iunit)


end subroutine open_gidmsh2_geometry
!
!
!
!
!
subroutine open_gidmsh2_geometry_mem(filename,ntri,npoints)
  implicit none

  !List of calling arguments
  character (len=*), intent(in) :: filename
  integer ( kind = 8 ), intent(out) :: npoints,ntri

    ! List of local variable
  integer ( kind = 8 ) :: io,iunit,i,nel,ind,ielem,ntag,ierror
  character ( len = 1000) :: cline

  iunit = 18
  open(UNIT=18, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ierror)
  read(iunit,*) ind,ind,ind,ind,ntri
  npoints = ntri*6
  close(iunit)

return
end subroutine open_gidmsh2_geometry_mem
!
!
!
!
!
