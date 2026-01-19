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
subroutine open_gmsh2_geometry(filename,npatches,norders,ixyzs, &
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
  integer *8 n,nomax,numnodes,numelem,ntri,ntag,npols,norder,nmax,nel
  integer *8 lll,m,l,j,iunit,istart,ipt,io,inode,info,i,ielem,ind
  real *8 x,y,z,rr
  character ( len = 1000) :: cline
 
!
!  Gather all node and vals2coefs (umats), and coefs2vals (pmat,dmat,vmat) 
!  matrices
!
 
  nomax = 4
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
        vmat(j,l) = pols(j)
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
        vmat(j,l) = pols(j)
      enddo
    enddo
    call dinverse(npols,vmat,info,umat)
    do j=1,npols
      do l=1,npols
        umatvrall(l,j,i) = umat(l,j)
      enddo
    enddo
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

  do
    read(iunit, *, iostat=io) cline

    if (io .ne. 0) exit
    !print *, 'Read line: ', trim(cline)

    if (trim(cline) .eq. '$Nodes') then
      print *, 'Reading nodes . . . '
      read(iunit,*) numnodes
      print *, 'Number of nodes = ', numnodes
      print *
      
      allocate(xyzs(3,numnodes))
      do i = 1,numnodes
        read (iunit,*) ind, x, y, z
        xyzs(1,i) = x
        xyzs(2,i) = y
        xyzs(3,i) = z
      end do

    end if

    if (trim(cline) .eq. '$Elements') then
      print *, 'Reading elements . . . '
      read(iunit,*) numelem
      print *, 'Number of elements = ', numelem
      print *

      ntri = 0
      do i = 1,numelem
        
        read(iunit, '(a)', iostat=io) cline
        read (cline,*) ind, ielem, ntag
        norder = -1
        if(ielem.eq.2) then
          ntri = ntri + 1
          norder = 1
        endif
        if(ielem.eq.9) then
          ntri = ntri + 1
          norder = 2
        endif
        if(ielem.eq.21) then
          ntri = ntri + 1
          norder = 3
        endif
        if(ielem.eq.23) then
          ntri = ntri + 1
          norder = 4
        endif

        if (norder .gt. 0) then

          ntri = ntri + 1
          npols = (norder+1)*(norder+2)/2
          nel = npols

          lll= 1+1+1+ntag+nel
          allocate(element(lll))
          read(cline,*) element
!
!  Temporarily store xyz coordinates for the triangle in 
!  srcvals (even though it is still gmsh nodes) srcvals 
!  will be overriden later
!
          do j = 1,npols
            inode = element(j+3+ntag)
            srcvals(1:3,ixyzs(ntri)+j-1) = xyzs(1:3,inode)
          end do

          istart = ixyzs(ntri)
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
                pmatall(m,l,norder)*srccoefs(1:3,istart+l-1)
              srcvals(4:6,istart+l-1) = srcvals(4:6,istart+l-1) + &
                dmatall(1,m,l,norder)*srccoefs(1:3,istart+l-1)
              srcvals(7:9,istart+l-1) = srcvals(7:9,istart+l-1) + &
                dmatall(2,m,l,norder)*srccoefs(1:3,istart+l-1)
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
          norders(ntri) = norder
          iptype(ntri) = 1
          ixyzs(ntri+1) = ixyzs(ntri) + npols


          deallocate(element)

        end if
      end do
    end if
  end do

  close(iunit)


end subroutine open_gmsh2_geometry
!
!
!
!
!
subroutine open_gmsh2_geometry_mem(filename,ntri,npoints)
  implicit none

  !List of calling arguments
  character (len=*), intent(in) :: filename
  integer ( kind = 8 ), intent(out) :: npoints,ntri

    ! List of local variable
  integer ( kind = 8 ) :: io,iunit,i,nel,ind,ielem,ntag
  integer *8 numelem,ierror
  character ( len = 1000) :: cline

  iunit = 18
  open(UNIT=18, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ierror)
  do
    read(iunit, *, iostat=io) cline

    if (io .ne. 0) exit
    !print *, 'Read line: ', trim(cline)

    if (trim(cline) .eq. '$Elements') then
      print *, 'Reading elements . . . '
      read(iunit,*) numelem
      print *, 'Number of elements = ', numelem
      print *

      ntri = 0
      do i = 1,numelem
        
        read(iunit, '(a)', iostat=io) cline
        read (cline,*) ind, ielem, ntag

        if(ielem.eq.2) then
          ntri = ntri + 1
          npoints = npoints + 3
        endif
        if(ielem.eq.9) then
          ntri = ntri + 1
          npoints = npoints + 6
        endif
        if(ielem.eq.21) then
          ntri = ntri + 1
          npoints = npoints + 10
        endif
        if(ielem.eq.23) then
          ntri = ntri + 1
          npoints = npoints + 15
        endif
      end do
    end if
  end do
  close(iunit)

return
end subroutine open_gmsh2_geometry_mem
!
!
!
!
!
subroutine load_uv_gmsh2_all(nomax,nmax,uv)
!
!  This subroutine loads the discretization nodes
!  for a gmsh discretization with order of discretization
!  = (1,2,3,or nomax <=4)
!
!  Input arguments:
!    - nomax: integer *8
!        max order of discretization nodes to be read (must be less than or 
!        equal to 4 currently)
!    - nmax: integer *8
!        max number of points to be read, must be = (nomax+1)*(nomax+2)/2
!    
!  Output arugments:
!    - uv: real *8 (2,nmax)
!        uv coordinates of the gmsh nodes on the standard triangles
!
!
  implicit none
  integer *8 nomax,nmax
  real *8 uv(2,nmax,*)


  uv(1,1,1) = 0
  uv(2,1,1) = 0
  
  uv(1,2,1) = 1
  uv(2,2,1) = 0

  uv(1,3,1) = 0
  uv(2,3,1) = 1

  if(nomax.eq.1) return

  uv(1,1,2) = 0
  uv(2,1,2) = 0
  
  uv(1,2,2) = 1
  uv(2,2,2) = 0

  uv(1,3,2) = 0
  uv(2,3,2) = 1

  uv(1,4,2) = 0.5d0
  uv(2,4,2) = 0

  uv(1,5,2) = 0.5d0
  uv(2,5,2) = 0.5d0

  uv(1,6,2) = 0.0d0
  uv(2,6,2) = 0.5d0

  if(nomax.eq.2) return
  

  uv(1,1,3) = 0
  uv(2,1,3) = 0
  
  uv(1,2,3) = 1
  uv(2,2,3) = 0

  uv(1,3,3) = 0
  uv(2,3,3) = 1

  uv(1,4,3) = 1.0d0/3.0d0
  uv(2,4,3) = 0

  uv(1,5,3) = 2.0d0/3.0d0
  uv(2,5,3) = 0

  uv(1,6,3) = 2.0d0/3.0d0
  uv(2,6,3) = 1.0d0/3.0d0

  uv(1,7,3) = 1.0d0/3.0d0
  uv(2,7,3) = 2.0d0/3.0d0

  uv(1,8,3) = 0
  uv(2,8,3) = 2.0d0/3.0d0

  uv(1,9,3) = 0
  uv(2,9,3) = 1.0d0/3.0d0

  uv(1,10,3) = 1.0d0/3.0d0
  uv(2,10,3) = 1.0d0/3.0d0

  if(nomax.eq.3) return

  uv(1,1,4) = 0
  uv(2,1,4) = 0
  
  uv(1,2,4) = 1
  uv(2,2,4) = 0

  uv(1,3,4) = 0
  uv(2,3,4) = 1

  uv(1,4,4) = 1.0d0/4.0d0
  uv(2,4,4) = 0

  uv(1,5,4) = 2.0d0/4.0d0
  uv(2,5,4) = 0

  uv(1,6,4) = 3.0d0/4.0d0
  uv(2,6,4) = 0 

  uv(1,7,4) = 3.0d0/4.0d0
  uv(2,7,4) = 1.0d0/4.0d0

  uv(1,8,4) = 2.0d0/4.0d0
  uv(2,8,4) = 2.0d0/4.0d0

  uv(1,9,4) = 1.0d0/4.0d0
  uv(2,9,4) = 3.0d0/4.0d0

  uv(1,10,4) = 0
  uv(2,10,4) = 3.0d0/4.0d0

  uv(1,11,4) = 0
  uv(2,11,4) = 2.0d0/4.0d0

  uv(1,12,4) = 0
  uv(2,12,4) = 1.0d0/4.0d0

  uv(1,13,4) = 1.0d0/4.0d0
  uv(2,13,4) = 1.0d0/4.0d0

  uv(1,14,4) = 2.0d0/4.0d0
  uv(2,14,4) = 1.0d0/4.0d0

  uv(1,15,4) = 1.0d0/4.0d0
  uv(2,15,4) = 2.0d0/4.0d0

end subroutine load_uv_gmsh2_all



subroutine load_uvw_vr_all(nomax,nmax,uv,w)
!
!  This subroutine loads the discretization nodes
!  for a vioreanu rokhlin discretization with max order of discretization
!  nomax. This routine also loads the quadrature weights
!
!  Input arguments:
!    - nomax: integer *8
!        max order of discretization nodes to be read
!    - nmax: integer *8
!        max number of points to be read, must be = (nomax+1)*(nomax+2)/2
!    
!  Output arugments:
!    - uv: real *8 (2,nmax,nomax)
!        uv coordinates of the vioreanu-rokhlin (vr) nodes on the standard 
!        triangles
!    - w: real *8 (nmax,nomax)
!        qaudrature weights for smooth functions for the vr nodes 
!        on the standard triangle for all orders <=nomax
!
!
  implicit none
  integer *8 nomax,nmax
  real *8 uv(2,nmax,nomax),w(nmax,nomax)
  integer *8 i,norder,npols,j

  do i=1,nomax
    do j=1,nmax
      uv(1,j,i) = 0
      uv(2,j,i) = 0
      w(j,i) = 0
    enddo
    npols = (i+1)*(i+2)/2
    call get_vioreanu_nodes_wts(i,npols,uv(1,1,i),w(1,i)) 
  enddo

end subroutine load_uvw_vr_all

