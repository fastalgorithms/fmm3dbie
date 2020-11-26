subroutine write_wtorus(fname,radii,scales,nosc,nu,nv,norder,pmax)
!
!f2py intent(in) fname,radii,scales,nosc,nu,nv,norder
!f2py intent(out) pmax
!
  implicit real *8 (a-h,o-z)
  character (len=*), intent(in) :: fname
  real *8, intent(in), target :: radii(3),scales(3)
  integer, intent(in) :: nu,nv,nosc,norder
  real *8, intent(out) :: pmax

  real *8, target :: p4(1)
  real *8, pointer :: ptr1,ptr2,ptr3,ptr4 
  real *8, allocatable, target :: triaskel(:,:,:)
  real *8, allocatable :: pdis(:),wts(:)
  real *8, allocatable :: srccoefs(:,:),srcvals(:,:)
  integer, allocatable :: norders(:),iptype(:),ixyzs(:)
  

  procedure (), pointer :: xtri_geometry

  external xtri_wtorus_eval

  real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts0(:)
 
  npols = (norder+1)*(norder+2)/2
  allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols), &
    wts0(npols))
  
  call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts0)

  npatches = nu*nv*2
  npts = npatches*npols
  done = 1
  pi = atan(done)*4
  umin = 0
  umax = 2*pi
  vmin = 0
  vmax = 2*pi
  
  allocate(triaskel(3,3,npatches))
  nover = 0
  call xtri_rectmesh_ani(umin,umax,vmin,vmax,nu,nv,nover,npatches, &
    npatches,triaskel)

  p4(1) = nosc + 0.0d0
  
  ptr1 => triaskel(1,1,1)
  ptr2 => radii(1)
  ptr3 => scales(1)
  ptr4 => p4(1)
  
  xtri_geometry => xtri_wtorus_eval

  allocate(srcvals(12,npts),srccoefs(9,npts))
  
  call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,npols, &
    uvs,umatr,srcvals,srccoefs)

  
  open(unit=33,file=trim(fname))
  write(33,*) norder
  write(33,*) npatches

  do l=1,12
    do i=1,npts
      write(33,*) srcvals(l,i) 
    enddo
  enddo

  close(33)

  allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))
  do i=1,npatches
    norders(i) = norder
    iptype(i) = 1
    ixyzs(i) = (i-1)*npols+1
  enddo
  ixyzs(npatches+1) = npts+1

  allocate(pdis(npatches),wts(npts))
  call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)
  do i=1,npatches
    pdis(i) = 0
  enddo

  pmax = 0.0d0

  call get_patch_distortion(npatches,norders,ixyzs,iptype,npts, &
   srccoefs,srcvals,wts,pdis)

  pmax = maxval(pdis)

end subroutine write_wtorus



subroutine get_wtorus_geom(radii,scales,nosc,nu,nv,npatches, &
  norder,npts,norders,ixyzs,iptype,srcvals,srccoefs,wts)
!
!f2py intent(in) radii,scales,nosc,nu,nv,npatches,norder,npts
!f2py intent(out) norders,ixyzs,iptype,srcvals,srccoefs,wts
!

  implicit real *8 (a-h,o-z)
  real *8, intent(in), target :: radii(3),scales(3)
  integer, intent(in) :: nu,nv,nosc,norder
  integer, intent(in) :: npatches,npts

  integer, intent(out) :: norders(npatches),ixyzs(npatches+1)
  integer, intent(out) :: iptype(npatches)

  real *8, intent(out) :: srcvals(12,npts),srccoefs(9,npts),wts(npts)

  real *8, target :: p4(1)
  real *8, pointer :: ptr1,ptr2,ptr3,ptr4 
  real *8, allocatable, target :: triaskel(:,:,:)
  
  external xtri_wtorus_eval

  procedure (), pointer :: xtri_geometry

  real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts0(:)
 
  npols = (norder+1)*(norder+2)/2
  allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols), &
    wts0(npols))
  
  call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts0)

  done = 1
  pi = atan(done)*4
  umin = 0
  umax = 2*pi
  vmin = 0
  vmax = 2*pi
  
  allocate(triaskel(3,3,npatches))
  nover = 0
  call xtri_rectmesh_ani(umin,umax,vmin,vmax,nu,nv,nover,npatches, &
    npatches,triaskel)

  p4(1) = nosc + 0.0d0
  
  ptr1 => triaskel(1,1,1)
  ptr2 => radii(1)
  ptr3 => scales(1)
  ptr4 => p4(1)
  
  xtri_geometry => xtri_wtorus_eval
  
  call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,npols, &
    uvs,umatr,srcvals,srccoefs)

  do i=1,npatches
    norders(i) = norder
    iptype(i) = 1
    ixyzs(i) = (i-1)*npols+1
  enddo
  ixyzs(npatches+1) = npts+1

  call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)

end subroutine get_wtorus_geom
!
!
!
!
!
subroutine get_sphere_geom(nref,npatches, &
  norder,npts,norders,ixyzs,iptype,srcvals,srccoefs,wts)
!
!f2py intent(in) radii,scales,nosc,nu,nv,npatches,norder,npts
!f2py intent(out) norders,ixyzs,iptype,srcvals,srccoefs,wts
!

  implicit real *8 (a-h,o-z)
  integer, intent(in) :: nref,norder
  integer, intent(in) :: npatches,npts

  integer, intent(out) :: norders(npatches),ixyzs(npatches+1)
  integer, intent(out) :: iptype(npatches)

  real *8, intent(out) :: srcvals(12,npts),srccoefs(9,npts),wts(npts)

  real *8, target :: p4(1),p2(1),p3(1)
  real *8, pointer :: ptr1,ptr2,ptr3,ptr4 
  real *8, allocatable, target :: triaskel(:,:,:)
  
  external xtri_wtorus_eval,xtri_sphere_eval

  procedure (), pointer :: xtri_geometry

  real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts0(:)
  integer, allocatable :: isides(:)
 
  npols = (norder+1)*(norder+2)/2
  allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols), &
    wts0(npols))
  
  call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts0)

  itype = 2
  allocate(triaskel(3,3,npatches))
  allocate(isides(npatches))
  npmax = npatches
  ntri = 0
  call xtri_platonic(itype, nref, npmax, ntri, triaskel, isides)

  xtri_geometry => xtri_sphere_eval
  ptr1 => triaskel(1,1,1)
  ptr2 => p2(1)
  ptr3 => p3(1)
  ptr4 => p4(1)


  
  call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,npols, &
    uvs,umatr,srcvals,srccoefs)

  do i=1,npatches
    norders(i) = norder
    iptype(i) = 1
    ixyzs(i) = (i-1)*npols+1
  enddo
  ixyzs(npatches+1) = npts+1

  call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)

end subroutine get_sphere_geom


