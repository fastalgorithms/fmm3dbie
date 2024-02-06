!    This file contains the following user callable subroutines
!      surf_fun_error - estimate the resolution of a function on
!         a surface
!
!      get_centroid_rads - compute the centroid and bounding sphere
!          radii for a collection of patches
!         
!      get_centroid_rads_guru - compute the centroid and bounding
!          sphere radius for a given triangular patch
!
!      oversample_geom - oversample xyz,dxyz/du,dxyz/dv on the surface
!
!      oversample_fun_surf - oversample functions defined on a surface
!
!      oversample_fun_guru - oversample function on a triangle
!
!      refine_patches_list - refine a given list of patches, assumes
!         that the input source arrays are long enough to accomodate
!         the additional patches to be added
!
!      surf_vals_to_coefs - given function values defined on a patch
!         convert to its basis expansion coeffs 
!
!      vals_to_coefs_guru - given function values defined on a patch
!         convert to its basis expansion coeffs 
!
!      get_qwts - compute smooth quadrature weights on the surface
!
!      get_near_far_split_pt - split the near field into the near-near
!         field (done via adaptive integration) and the near-far field
!         (done via oversampled quadrature (more oversampling than the
!         far-far field))
!
!      get_patch_id_uvs - for all boundary points, get patch id they are on
!          and u,v location on patch
!
!      get_patch_distortion - estimate shape distortion parameter for a given
!          collection of patches
!
!      get_nximat - estimate number of elements in the union of 
!        interpolation matrices for all patches
!
!      get_ximats - get all interpolation matrices for the surface
!   
!      get_first_fundamental_form - get the first fundamental form at all
!        points
!
!      get_inv_first_fundamental_form - get inverse of first fundamental 
!         form
!
!      get_surf_grad - compute surface gradient of a scalar function
!
!      get_surf_grad_fast - compute surface gradient of a scalar function 
!        (With precomputed inverse of first fundamental form)
!
!      get_surf_uv_grad - compute the uv gradient of a function defined
!         on a surface, i.e. \partial_{u}f, \partial_{v}f
!
!      get_surf_div - compute the surface divergence of a vector field
!         defined on a surface
!
!      col_ind_to_patch_node_ind - convert a list of column
!         indices to a list of patch and node indices
!
!      get_second_fundamental from - get the second fundamental form
!         at the discretization nodes on the surface
!
!      get_mean_curvature - compute the mean curvature at all discretization
!         points on the surface
!      
!
!         

subroutine surf_fun_error(nd,npatches,norders,ixyzs,iptype, &
   npts,dvals,wts,errp,errm)

!
!  This subroutine computes the error in resolving
!  a function of the surface whose patchwise
!  basis expansion coefficients are provided
!
!  nd - number of functions
!  npatches - number of patches
!  norders(npatches) - order of discretization
!  ixyzs(npatches+1) - location where expansion coeffs for patch i start
!  iptype(npatches) - type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!  dvals(nd,npts) - values of the function at discretization
!                   nodes
!  wts(npts) - quadrature weights for integrating smooth functions on
!              the surface
!  
!
!  output:
!  errp(npatches) - error per patch
!  errm - max of errp
!

  implicit real *8 (a-h,o-z)
  integer norders(npatches),ixyzs(npatches+1)
  integer iptype(npatches)
  real *8 dvals(nd,npts),errp(npatches),errm
  real *8 wts(npts)
  real *8, allocatable :: dcoefs(:,:),dtail(:,:),uvs(:,:),pols(:)
  integer, allocatable :: itailcoefs(:)
  real *8 rl2s,rl2tails

  allocate(dcoefs(nd,npts),dtail(nd,npts))

  do i=1,npts
    do idim=1,nd
      dcoefs(idim,i) = 0
      dtail(idim,i) = 0
    enddo
  enddo

  call surf_vals_to_coefs(nd,npatches,norders,ixyzs,iptype,npts, &
   dvals,dcoefs)
  
  
  errm = 0
  rl2s = 0
  rl2tails = 0


  nmax = 20
  npmax = (nmax+1)*(nmax+2)

  allocate(uvs(2,npmax),pols(npmax),itailcoefs(npmax))
  do ip=1,npatches
    norder = norders(ip)
    istart0 = ixyzs(ip)
    npols = ixyzs(ip+1)-ixyzs(ip)
    call get_disc_nodes(norder,npols,iptype(ip),uvs)
    call get_tail_coefs(norder,npols,iptype(ip),itailcoefs,ntailcoefs)



    if(iptype(ip).eq.1) then
      do i=1,npols
        ipt = ixyzs(ip)+i-1
        call get_basis_pols(uvs(1,i),norder,npols,iptype(ip),pols)
        do idim=1,nd
          dtail(idim,ipt) = 0
          do jj=1,ntailcoefs
            j = itailcoefs(jj)
            ic = ixyzs(ip)+j-1
            dtail(idim,ipt) = dtail(idim,ipt) + pols(j)*dcoefs(idim,ic)
          enddo
        enddo
      enddo
    endif

    errp(ip) = 0
    do idim=1,nd
      do i=1,npols
        ipt = ixyzs(ip)+i-1
        errp(ip) = errp(ip) + dtail(idim,ipt)**2*wts(ipt)
        rl2s = rl2s + dvals(idim,ipt)**2*wts(ipt)
        rl2tails = rl2tails + dtail(idim,ipt)**2*wts(ipt)
      enddo
    enddo
    errp(ip) = sqrt(errp(ip))
    if(errp(ip).gt.errm) errm = errp(ip)
  enddo
  
  rl2s = sqrt(rl2s)
  rl2tails = sqrt(rl2tails)

  do ip=1,npatches
    errp(ip) = errp(ip)/errm*rl2tails/rl2s
  enddo

  errm = rl2tails/rl2s


end subroutine surf_fun_error
!
!
!
!
!




subroutine get_centroid_rads(npatches,norders,ixyzs,iptype,npts, &
    srccoefs,cms,rads)
!
!   this subroutine computes the centroid of each patch and the radius 
!   of the bounding sphere centered at the centroid
!
!   Input:
!      npatches - integer
!         number of patches
!      norders - integer
!         order of discretization of each patches
!      ixyzs - integer(npatches+1)
!            ixyzs(i) denotes the starting location in srccoefs,
!             and srcvals array corresponding to patch i
!   
!      iptype - integer(npatches)
!         type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!
!     npts - integer
!         total number of discretization points
!      srccoefs - real *8 (9,npts)
!         basis expansions of the geometry info
!           x,y,z,,dx/du,dy/du,dz/du,dx/dv,dy/dv,dz/dv
!
!

  implicit none
  integer npatches,norders(npatches),npols,npts,ixyzs(npatches+1)
  integer iptype(npatches)
  real *8 srccoefs(9,npts),cms(3,npatches),rads(npatches)

  integer i,istart


  do i=1,npatches
    istart = ixyzs(i)
    npols = ixyzs(i+1)-ixyzs(i)
    call get_centroid_rads_guru(norders(i),npols, iptype(i), &
      srccoefs(1,istart),cms(1,i),rads(i))
  enddo
end subroutine get_centroid_rads
!
!
!
!

subroutine get_centroid_rads_guru(norder,npols,iptype,srccoefs,cms,rads) 
  implicit none
  integer norder,npols,iptype
  real *8 srccoefs(9,npols),cms(3),rads
  real *8 uv(2,4)
  real *8 xyz(3,4),rad
  real *8, allocatable :: pols(:,:)

  integer i,j,l,m,lpt,np,nv

  call get_boundary_vertices(iptype,uv,nv)
  

  allocate(pols(npols,nv))
  do i=1,nv
    call get_basis_pols(uv(1,i),norder,npols,iptype,pols(1,i))
  enddo

  do l=1,3
    cms(l) = 0
  enddo
  rads = 0

!
!     loop over all three vertices
!
  do j=1,nv
    do m=1,3
      xyz(m,j) = 0
    enddo
    
    do l=1,npols
      do m=1,3
        xyz(m,j) = xyz(m,j) + srccoefs(m,l)*pols(l,j)
      enddo
    enddo
  enddo

  do m=1,3
    do l=1,nv
      cms(m) = cms(m) + xyz(m,l)
    enddo
    cms(m) = cms(m)/nv
  enddo
!
!    compute radius of bounding sphere
!
  do m=1,nv
    rad = (cms(1)-xyz(1,m))**2 + (cms(2)-xyz(2,m))**2 + &
        (cms(3)-xyz(3,m))**2
    if(rad.ge.rads) rads = rad
  enddo
  rads = sqrt(rads)

  return
end subroutine get_centroid_rads_guru
!
!
!
!
subroutine oversample_geom(npatches,norders,ixyzs,iptype, &
    npts,srccoefs,srcvals,nfars,ixyzso,nptso,srcover)
!
!
!  This subroutine oversamples geometry information
!  given expansion coeffs of geometry info on patches
!
!
!  Extreme care must be taken when an oversampled node is identical
!  to a discretization node. The routine ensures that in this event
!  the source information is copied over from the original
!  discretization
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        discretization order of patches
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!    - nfars: integer(npatches)
!        oversampled order of patch i
!    - ixyzso: integer(npatches+1)
!        starting location of oversampled points on patch i
!    - nptso: integer
!        total number of oversampled points
!
!  Output arguments:   
!    
!    - srcover: double precision(12,nptso)
!        Oversampled geometry information
!
!----------------------------------------------
!

  implicit none
  integer, intent(in) :: npatches,norders(npatches),ixyzs(npatches+1)
  integer, intent(in) :: iptype(npatches),ixyzso(npatches+1)
  integer, intent(in) :: nfars(npatches)
  integer, intent(in) :: npts,nptso
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
  real *8, intent(out) :: srcover(12,nptso)
  integer nfar,norder
  real *8 tmp(3),rr
  

  real *8, allocatable :: pmat(:,:),pols(:)
  real *8, allocatable :: uvs(:,:),umat(:,:),vmat(:,:),wts(:)
  character *1 transa,transb
  real *8 alpha, beta

  integer i,istart,istartover,nfar_pols,j,jpt,npols,l
  integer n1,n2,ipt,i1,i2
  integer ic,nclash,nmax,npmax
  integer, allocatable :: iclash(:),iclashfar(:)


  transa = 'n'
  transb = 'n'
  alpha = 1
  beta = 0

  nmax = maxval(norders)
  npmax = (nmax+1)*(nmax+1)
  allocate(iclash(npmax),iclashfar(npmax))
  

  do i=1,npatches
    nfar = nfars(i)
    norder = norders(i)

    if(nfar.ne.norder) then
      nfar_pols = ixyzso(i+1)-ixyzso(i)
      npols = ixyzs(i+1)-ixyzs(i)
      allocate(uvs(2,nfar_pols))
      call get_disc_nodes(nfar,nfar_pols,iptype(i),uvs)


      allocate(pmat(npols,nfar_pols))
    
      do j=1,nfar_pols
        call get_basis_pols(uvs(1,j),norder,npols,iptype(i),pmat(1,j))
      enddo

      istart = ixyzs(i)
      istartover = ixyzso(i)
      call dgemm_guru(transa,transb,9,nfar_pols,npols,alpha, &
          srccoefs(1,istart),9,pmat,npols,beta,srcover(1,istartover),12)
      do j=1,nfar_pols
        jpt = istartover+j-1
        call cross_prod3d(srcover(4,jpt),srcover(7,jpt),tmp)
        rr = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
        srcover(10,jpt) = tmp(1)/rr
        srcover(11,jpt) = tmp(2)/rr
        srcover(12,jpt) = tmp(3)/rr
      enddo
      
      nclash = 0
      call get_clashing_indices(norder,npols,iptype(i),nfar,nfar_pols, &
        iclash,iclashfar,nclash)

!
!    Fix clash in coordinates across orders
!
      do ic=1,nclash
        ipt = ixyzs(i)+iclash(ic)-1
        jpt = ixyzso(i) + iclashfar(ic)-1
        do l=1,12
          srcover(l,jpt) = srcvals(l,ipt)
        enddo
      enddo
      deallocate(uvs,pmat)
    else
      npols = ixyzs(i+1)-ixyzs(i)
      istart = ixyzs(i)
      istartover = ixyzso(i)
      do j=1,npols
        do l=1,12
          srcover(l,j+istartover-1) = srcvals(l,j+istart-1)
        enddo
      enddo
    endif
  enddo
end subroutine oversample_geom





subroutine oversample_fun_surf(nd,npatches,norders,ixyzs,iptype,npts, &
   u,nfars,ixyzso,nptso,uover)
!
!  This subroutine oversamples a function defined on a surface
!
!  input
!    nd - integer
!      number of functions
!    npatches - integer
!      number of patches on the surface
!    norders - integer(npatches)
!      order of discretization of patch i
!    ixyzs - integer(npatches+1)
!      starting location of points on patch i
!    iptype - integer(npatches)
!      type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!    npts - integer
!      total number of points on the surface
!    u - real *8 (nd,npts)
!      functions at the discretization points on the surface
!    nfars - integer(npatches)
!      oversampled order of patch i
!    ixyzso - integer(npatches+1)
!      starting location of oversampled points on patch i
!    nptso - integer
!      total number of oversampled points
!  
!
!  output
!    uover - real *8 (nd,nptso)
!       oversampled function
!    
   

  implicit none
  integer nd,npatches,norders(npatches),ixyzs(npatches+1)
  integer iptype(npatches),npts
  real *8 u(nd,npts)
  integer nfars(npatches),ixyzso(npatches+1),nptso
  real *8 uover(nd,nptso)
  integer i,istart,istarto,npols,npolso

  do i=1,npatches
    
    istart = ixyzs(i)
    istarto = ixyzso(i)
    npols = ixyzs(i+1)-ixyzs(i)
    npolso = ixyzso(i+1) - ixyzso(i)
    call oversample_fun_guru(nd,norders(i),npols,iptype(i),u(1,istart),&
      nfars(i),npolso,uover(1,istarto))
  enddo

end subroutine oversample_fun_surf
!
!
!
!
subroutine oversample_fun_guru(nd,norder,npols,iptype,u,nfar,nfar_pols, &
    uover)
!
!  This subroutine oversample a function on the standard triangle
!
!  input
!    nd - integer
!       number of functions
!    norder - integer
!       order of original discretization
!    npols - integer
!       number of discretization points corresponding to norder
!    iptype - integer
!       type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!    u - real *8 (nd,npols)
!       function values tabulated at original grid
!    nfar - integer
!       oversampled order of discretization
!    nfar_pols - integer
!       number of discretization nodes corresponding to oversampled
!         order
!    
!  output
!    uover - real *8 (nd,nfar_pols)
!      oversampled function values
!


  implicit none
  integer norder,npover,nd,npols,nfar
  real *8 u(nd,npols),uover(nd,nfar_pols)
  real *8, allocatable :: ucoefs(:,:)
  real *8, allocatable :: pmat(:,:),pols(:)
  real *8, allocatable :: uvs(:,:),umat(:,:),vmat(:,:),wts(:)
  real *8, allocatable :: umat0(:,:),uvs0(:,:),vmat0(:,:),wts0(:)
  character *1 transa,transb
  integer iptype
  real *8 alpha, beta
  integer i,nfar_pols


  allocate(uvs(2,nfar_pols))
  call get_disc_nodes(nfar,nfar_pols,iptype,uvs)

  allocate(ucoefs(nd,npols))

  allocate(umat0(npols,npols),vmat0(npols,npols),uvs0(2,npols))
  allocate(wts0(npols))


  call get_disc_exps(norder,npols,iptype,uvs0,umat0,vmat0,wts0)


!
!  compute ucoefs
!

  transa = 'n'
  transb = 't'
  alpha = 1
  beta = 0
  call dgemm_guru(transa,transb,nd,npols,npols,alpha,u,nd,umat0,npols,&
     beta,ucoefs,nd)

  allocate(pmat(npols,nfar_pols))
    
  do i=1,nfar_pols
    call get_basis_pols(uvs(1,i),norder,npols,iptype,pmat(1,i))
  enddo
    
  transa = 'n'
  transb = 'n'
  call dgemm_guru(transa,transb,nd,nfar_pols,npols,alpha, &
        ucoefs,nd,pmat,npols,beta,uover,nd)
 
end subroutine oversample_fun_guru
!
!
!
!
!
subroutine refine_patches_list(npatches, npatmax, norders, ixyzs, &
   iptype, npts, nptmax, srccoefs, srcvals, nlist, ilist, &
   npatches_out, npts_out, ier)
!
!
!  This subroutine refines a given set of patches on a discretized 
!  surface. The routine assumes that all the arrays are long 
!  enough to accomodate the refined surface, i.e. npatmax must
!  be at least npatches + 3*nlist, and nptmax must be atleast
!
!  nptmax >= npts + \sum_{i} (ixyzs(ilist(i)+1) - ixyzs(ilist(i)))*3
!
!  If not, the routine will return without refining any patch, and
!  an error code of 4.
!
!  Developer notes: 
!  When refining, the first child patch in the tree convention used
!  in adaptive integration routines will be left in place of the 
!  original patch, and the other 3 patches will be appended to
!  the end of the list
!
!  Input arguments:
!    - npatches: integer
!        number of patches describing the discretized surface
!    - npatmax: integer
!        maximum number of patches and size of norders, ixyzs (-1),
!        iptype on input
!    - norders: integer(npatmax)
!        norders(1:npatches) is the discretization order of patches
!    - ixyzs: integer(npatmax+1)
!        ixyzs(1:npatches+1) is starting location of points on patch i
!    - iptype: integer(npatmax)
!        iptype(1:npatches) type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!    - npts: integer
!        total number of points describing the discretized surface
!    - nptmax: integer
!        maximum number of points, and size of second dimension of
!        srccoefs, and srcvals arrays
!    - srccoefs: double precision (9,nptmax)
!        srccoefs(1:9,1:npts) are the koornwinder expansion coefs of 
!        geometry info
!    - srcvals: double precision (12,nptmax)
!        srcvals(1:12,1:npts) are xyz, dxyz/du,dxyz/dv, normals at 
!        all nodes
!    - nlist: integer
!        number of patches of to be refined
!    - ilist: integer(nlist)
!        list of patches to be refined
!
!  Output arguments
!    - npatches_out: integer
!        number of patches on the refined surface 
!    - norders: integer(npatmax)
!        norders(1:npatches_out) is the discretization order of patches
!    - ixyzs: integer(npatmax+1)
!        ixyzs(1:npatches_out+1) is starting location of points on patch i
!    - iptype: integer(npatmax)
!        iptype(1:npatches_out) type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!    - npts_out: integer
!        total number of points describing the refined surface
!    - srccoefs: double precision (9,nptmax)
!        srccoefs(1:9,1:npts_out) are the koornwinder expansion coefs 
!        of geometry info
!    - srcvals: double precision (12,nptmax)
!        srcvals(1:12,1:npts_out) are xyz, dxyz/du,dxyz/dv, normals at 
!        all nodes
!    - ier: integer
!        error code, 
!        * ier = 0, implies successful refinement,
!        * ier = 4, implies that either nptmax or npatmax is too small
!     
  implicit none
  integer, intent(in) :: npatches, npatmax, npts, nptmax
  integer, intent(inout) :: norders(npatmax), ixyzs(npatmax+1)
  integer, intent(inout) :: iptype(npatmax)
  real *8, intent(inout) :: srccoefs(9,nptmax), srcvals(12,nptmax)
  integer, intent(in) :: nlist, ilist(nlist)
  integer, intent(out) :: npatches_out, npts_out, ier

!   temporary variables
  integer, allocatable :: ipat_starts(:), ipt_starts(:), npol_ilist(:)
  integer ntest, i, i1, ipat_start, ipt_start, j, k, l, m
  real *8, allocatable :: uvs_refine(:,:), umatr(:,:), uvs_base(:,:)
  real *8, allocatable :: vmatr(:,:), wts(:), pols(:)
  integer nordermax, npolmax, ipat, npols, jj, ll, ipats(4), jpat, kk
  integer iptype0, norder0, istart
  real *8 v0(2,3), v1(2,3), v2(2,3), v3(2,3), v4(2,3), tmp(3), rr


!   Test if arrays are large enough to handle the refined geometry
  if(npatmax.le.(npatches + 3*nlist)) then
    ier = 4
    npatches_out = npatches
    npts_out = npts
    return
  endif

   allocate(npol_ilist(nlist), ipat_starts(nlist), ipt_starts(nlist))
!$OMP PARALLEL DO DEFAULT(SHARED)
   do i=1,nlist
     npol_ilist(i) = 3*(ixyzs(ilist(i)+1) - ixyzs(ilist(i)))
     ipat_starts(i) = npatches + 3*(i-1) + 1
   enddo
!$OMP END PARALLEL DO

   call cumsum(nlist, npol_ilist, ipt_starts)
   ntest = npts + ipt_starts(nlist)
   
   i1 = ipt_starts(1)
!$OMP PARALLEL DO DEFAULT(SHARED)
   do i=1,nlist
     ipt_starts(i) = ipt_starts(i) - i1 + npts + 1
   enddo
!$OMP END PARALLEL DO

  
  if(nptmax.le.ntest) then
    ier = 4
    npatches_out = npatches
    npts_out = npts
    return
  endif

  ier = 0
  print *, ier
  call prinf('ipts_start=*',ipt_starts,nlist)
  call prinf('ntest=*',ntest,1)
  call prinf('ipat_starts=*',ipat_starts,nlist)

  npatches_out = npatches + 3*nlist
  npts_out = ntest


!  Start refining the geometry now  
!   

  nordermax = maxval(norders(1:npatches))
  npolmax = (nordermax+1)*(nordermax)
  allocate(uvs_base(2,npolmax))
  allocate(uvs_refine(2,4*npolmax), umatr(npolmax, npolmax))
  allocate(vmatr(npolmax, npolmax), wts(npolmax))
  allocate(pols(npolmax))

  

!
!  TODO: find unique order + patch type combo to precompute
!  required umatr's and then reuse them across
!  refinement work being done
!

!
!  fill out ixyzs,norders, iptype array indepedently 
!  to make sure all entries are in
!
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ipat,ipats,iptype0) &
!$OMP PRIVATE(norder0,npols)
   do i=1,nlist
     ipat = ilist(i)
     call get_npols(iptype(ipat), norders(ipat), npols)
     ipats(1) = ipat
     ipats(2) = ipat_starts(i)
     ipats(3) = ipat_starts(i)+1
     ipats(4) = ipat_starts(i)+2
     

     iptype0 = iptype(ipat)

     iptype(ipats(2)) = iptype0 
     iptype(ipats(3)) = iptype0 
     iptype(ipats(4)) = iptype0 

     norder0 = norders(ipat)
     norders(ipats(2)) = norder0 
     norders(ipats(3)) = norder0 
     norders(ipats(4)) = norder0 


! 
!  Note the extra shift by 1, since ixyzs(npatches+1)
!  was already filled out to be the correct start point
!
     ixyzs(ipats(2)+1) = ipt_starts(i) + npols
     ixyzs(ipats(3)+1) = ipt_starts(i) + 2*npols
     ixyzs(ipats(4)+1) = ipt_starts(i) + 3*npols
   enddo
!$OMP END PARALLEL DO   

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(uvs_refine, umatr) &
!$OMP PRIVATE(uvs_base, v0, v1, v2, v3, v4, ipat, npols, j) &
!$OMP PRIVATE(vmatr, wts, pols, istart, jj, ll, tmp, ipats) &
!$OMP PRIVATE(jpat, k, kk, l, m, rr)
  do i=1,nlist

!  Determine discretization nodes on refined patch  
     ipat = ilist(i)

     if(iptype(ipat).eq.1) then
       v0(1,1) = 0
       v0(2,1) = 0

       v0(1,2) = 1
       v0(2,2) = 0

       v0(1,3) = 0
       v0(2,3) = 1

       call gettrichildren(v0, v1, v2, v3, v4)
     elseif(iptype(ipat).eq.11.or.iptype(ipat).eq.12) then
       v0(1,1) = -1
       v0(2,1) = -1

       v0(1,2) = 1
       v0(2,2) = -1

       v0(1,3) = -1
       v0(2,3) = 1

       call getquadchildren(v0, v1, v2, v3, v4)
     endif
     call get_npols(iptype(ipat), norders(ipat), npols)
     call get_disc_exps(norders(ipat), npols, iptype(ipat), uvs_base, &
        umatr(1:npols,1:npols), vmatr, wts)

     if(iptype(ipat).eq.1) then
        call mapuv_tri(v1, npols, uvs_base, uvs_refine(1,1))
        call mapuv_tri(v2, npols, uvs_base, uvs_refine(1,1+npols))
        call mapuv_tri(v3, npols, uvs_base, uvs_refine(1,1+2*npols))
        call mapuv_tri(v4, npols, uvs_base, uvs_refine(1,1+3*npols))
     else if(iptype(ipat).eq.11.or.iptype(ipat).eq.12) then
        call mapuv_quad(v1, npols, uvs_base, uvs_refine(1,1))
        call mapuv_quad(v2, npols, uvs_base, uvs_refine(1,1+npols))
        call mapuv_quad(v3, npols, uvs_base, uvs_refine(1,1+2*npols))
        call mapuv_quad(v4, npols, uvs_base, uvs_refine(1,1+3*npols))
     endif

!
!   Fix iptype, ixyzs, norders
!
     ipats(1) = ipat
     ipats(2) = ipat_starts(i)
     ipats(3) = ipat_starts(i)+1
     ipats(4) = ipat_starts(i)+2




!
!   Fix srcvals
!
!
     do j=1,4*npols
       call get_basis_pols(uvs_refine(1,j), norders(ipat), npols, &
          iptype(ipat), pols)
       istart = ixyzs(ipat) 

       if(j.le.npols) then
         jj = istart+j-1
       else
         jj = ipt_starts(i) + j-npols-1
       endif

       do m=1,9
         srcvals(m, jj) = 0
       enddo

       do l=1,npols
          ll = istart+l-1
          do m=1,9
            srcvals(m,jj) = srcvals(m,jj) + srccoefs(m,ll)*pols(l)
          enddo
       enddo
       call cross_prod3d(srcvals(4,jj),srcvals(7,jj),tmp)
       rr = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
       srcvals(10,jj) = tmp(1)/rr
       srcvals(11,jj) = tmp(2)/rr
       srcvals(12,jj) = tmp(3)/rr
     enddo
!
!
!  Fix srccoefs
!
    do j=1,4
      jpat = ipats(j)
      do l=1,npols
        ll = ixyzs(jpat) + l-1
        do m=1,9
          srccoefs(m,ll) = 0.0d0 
        enddo

        do k=1,npols
          kk = ixyzs(jpat) + k-1
          do m=1,9
             srccoefs(m,ll) = srccoefs(m,ll) + umatr(l,k)*srcvals(m,kk) 
          enddo
        enddo
      enddo
    enddo

  enddo
!$OMP END PARALLEL DO

  return
end subroutine refine_patches_list


subroutine surf_vals_to_coefs(nd,npatches,norders,ixyzs,iptype,npts, &
   u,ucoefs)

!
!f2py intent(in) nd,npatches,norders,ixyzs,iptype,npts,u
!f2py intent(out) ucoefs
!
!
!  This subroutine oversamples a function defined on a surface
!
!  input
!    nd - integer
!      number of functions
!    npatches - integer
!      number of patches on the surface
!    norders - integer(npatches)
!      order of discretization of patch i
!    ixyzs - integer(npatches+1)
!      starting location of points on patch i
!    iptype - integer(npatches)
!      type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!    npts - integer
!      total number of points on the surface
!    u - real *8 (nd,npts)
!      functions at the discretization points on the surface
!
!  output
!    ucoefs - real *8 (nd,npts)
!       coeffs of funs in poly basis
!    
   

  implicit none
  integer nd,npatches,norders(npatches),ixyzs(npatches+1)
  integer iptype(npatches),npts
  real *8 u(nd,npts),ucoefs(nd,npts)
  integer i,istart,npols

  do i=1,npatches
    istart = ixyzs(i)
    npols = ixyzs(i+1)-ixyzs(i)
    call vals_to_coefs_guru(nd,norders(i),npols,iptype(i),u(1,istart),&
      ucoefs(1,istart))
  enddo

end subroutine surf_vals_to_coefs

!
!
subroutine vals_to_coefs_guru(nd,norder,npols,iptype,u,ucoefs)
!
!  This subroutine converts function to koornwinder expansion
!  coefs
!
!  input
!    nd - integer
!       number of functions
!    norder - integer
!       order of original discretization
!    npols - integer
!       number of discretization points corresponding to norder
!    iptype - integer
!       type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!    u - real *8 (nd,npols)
!       function values tabulated at original grid
!    
!  output
!    ucoefs - real *8 (nd,npols)
!       koornwinder expansion coefs
!


  implicit none
  integer norder,nd,npols,iptype
  real *8 u(nd,npols),ucoefs(nd,npols)
  real *8, allocatable :: umat0(:,:),uvs0(:,:),vmat0(:,:),wts0(:)
  character *1 transa,transb
  real *8 alpha, beta
  integer i



  allocate(umat0(npols,npols),vmat0(npols,npols),uvs0(2,npols))
  allocate(wts0(npols))


  call get_disc_exps(norder,npols,iptype,uvs0,umat0,vmat0,wts0)


!
!  compute ucoefs
!

  transa = 'n'
  transb = 't'
  alpha = 1
  beta = 0
  call dgemm_guru(transa,transb,nd,npols,npols,alpha,u,nd,umat0,npols,&
     beta,ucoefs,nd)
 
end subroutine vals_to_coefs_guru 
!
!
!
!
!
!
subroutine get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,qwts)

!
!f2py intent(in) npatches,norders,ixyzs,iptype,npts,srcvals
!f2py intent(out) qwts
!
!  this subroutine returns the quadrature weights
!  for integrating smooth functions on the surface
!
!
!  The surface is assumed to be discretizaed using
!  a collection of patches with the same order discretization
!  nodes
!
!  Input:
!    npatches - number of patches
!    norders(npatches) - order of discretization
!    ixyzs(npatches+1) - starting location of points on patch i
!    iptype(npatches) - type of patch
!    npts - total number of points on the surface
!    srcvals(12,npts) - surface info at the discretization nodes 
!
!  Output:
!    qwts - smooth quadrature weights at the discretization nodes
!
  implicit none
  integer norders(npatches),npatches,npts,npols,i,j,ipt
  integer ixyzs(npatches+1),iptype(npatches)
  real *8 srcvals(12,npts),qwts(npts),tmp(3)
  integer istart
  real *8, allocatable :: wts0(:)
  integer nmax,npmax
  
  nmax = maxval(norders(1:npatches))
  npmax = (nmax+1)*(nmax+1)
  allocate(wts0(npmax))

  do i=1,npatches
    istart = ixyzs(i)
    npols = ixyzs(i+1)-ixyzs(i)
    call get_disc_wts(norders(i),npols,iptype(i),wts0)
    do j=1,npols
      ipt = istart + j-1
      call cross_prod3d(srcvals(4,ipt),srcvals(7,ipt),tmp)
      qwts(ipt) = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)*wts0(j)
    enddo
  enddo
  deallocate(wts0)
end subroutine get_qwts

!
!
!
!
!
subroutine get_near_far_split_pt(nd,n,xyz,r,xyz0,nf,nn,xyz_f,xyz_n, &
    iind_f,iind_n)
!
!
!   This subroutine splits a given set of targets into
!   near points and far points with respect to a reference
!   point.
!
!   A target is labelled near if |t-xyz0|<=r \,, 
!   where xyz0 is the reference point, t is the target
!   
!
!   input
!     nd - integer
!        dimension of target array (must at least be 3, with the
!        first three components being the coordinates)
!     n - integer 
!       number of targets
!     xyz - real *8 (nd,n)
!       location of targets
!     r - real *8 
!       reference radius
!     xyz0 - real *8(3)
!        reference point
!     
!     output
!      nf - integer
!        number of far targets
!      nn - integer
!        number of near targets
!      xyz_f - real *8 (nd,n)
!        Collection of far targets. Note that
!        only the first (nd,nf) block will
!        be filled out
!      xyz_n - real *8 (nd,n)
!        Collection of near targets. Note that
!        only the first (nd,nn) block will
!        be filled out
!      iind_f - integer (n)
!         iind_f(i) denotes the location in xyz array
!         of the ith target in xyz_f.
!         Note that the first nf entries are filled out
!      iind_n - integer (n)
!         iind_n(i) denotes the location in xyz array
!         of the ith target in xyz_n.
!         Note that the first nn entries are filled out
!      
  implicit none
  integer nd
  integer n,nf,nn,iind_f(n),iind_n(n)
  real *8 xyz(nd,n),xyz0(3),r,rtmp,r2
  real *8 xyz_f(nd,n),xyz_n(nd,n)

  integer i,j

  r2 = r**2
  nf = 0
  nn = 0
  do i=1,n
    rtmp = (xyz(1,i)-xyz0(1))**2 + (xyz(2,i)-xyz0(2))**2 + &
       (xyz(3,i)-xyz0(3))**2
    if(rtmp.le.r2) then
      nn = nn +1

      do j=1,nd
        xyz_n(j,nn) = xyz(j,i)
      enddo

      iind_n(nn) = i
    else
      nf = nf +1
      do j=1,nd
        xyz_f(j,nf) = xyz(j,i)
      enddo
      iind_f(nf) = i
    endif
  enddo

end subroutine get_near_far_split_pt




subroutine get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts, &
    ipatch_id,uvs_pts)
!
!f2py intent(in) npatches,norders,ixyzs,iptype,npts
!f2py intent(out) ipatch_id,uvs_pts
!
  implicit none
  integer npatches,norders(npatches),ixyzs(npatches+1),iptype(npatches)
  integer npts,ipatch_id(npts)
  real *8 uvs_pts(2,npts)
  
  integer i,ip,j,npols


  do ip=1,npatches
    do j=ixyzs(ip),ixyzs(ip+1)-1
      ipatch_id(j) = ip
    enddo

    npols = ixyzs(ip+1)-ixyzs(ip)
    call get_disc_nodes(norders(ip),npols,iptype(ip),uvs_pts(1,ixyzs(ip)))
  enddo


end subroutine get_patch_id_uvs




subroutine get_patch_distortion(npatches,norders,ixyzs,iptype,npts,&
    srccoefs,srcvals,qwts,pdis)

!
!f2py intent(in) npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals
!f2py intent(in) qwts
!f2py intent(out) pdis
!
!
!   This subroutine estimates the patch distortion for a collection of 
!   patches:
!
!   For the right triangle, the patch distortion is defined as follows:
!
!   Let J denote the matrix whose columns are dxyz/du, and dxyz/dv
!
!   Let s1, s2 denote the largest and smallest singular values of the matrix
!     J 
!
!   Then the patch distortion is given by 
!      sqrt(\int_{T0} (s1/s2)^2 |det J| du dv/ \int_{T0} |det J| du dv)
!
!   
!   input:
!     npatches - number patches
!     norders(npatches) - order of discretization
!     ixyzs(npatches+1) - starting location of nodes on patch i
!     iptype(npatches) - type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!     npts - total number of discretization points
!     srcvals(12,npts) - xyz, dxyz/du,dxyz/dv, normals at all nodes
!     srccoefs(9,npts) - koornwinder expansion coeffs
!                         of geometry info
!     qwts(npts) - quadrature weights at discretization nodes
!
!     output:
!      pdis(npatches) - l2 distortion of each patch
!    
!
    
  implicit none
  integer npatches,norders(npatches),ixyzs(npatches+1),iptype(npatches)
  integer npts
  real *8 srccoefs(9,npts),srcvals(12,npts),qwts(npts),pdis(npatches)
  real *8 rtmp,srat,vtmp(3)
  integer i,istart,j,ipt,npols

  
  
  do i=1,npatches
    istart = ixyzs(i)
    npols = ixyzs(i+1)-ixyzs(i)
    pdis(i) = 0
    rtmp = 0
    do j=1,npols
      ipt = istart+j-1
      call get_local_asprat(srcvals(4,ipt),srcvals(7,ipt),srat)
      rtmp = rtmp + qwts(ipt)
      pdis(i) = pdis(i) + srat*qwts(ipt)
    enddo
    pdis(i) = sqrt(pdis(i)/rtmp)
  enddo
  return
end subroutine get_patch_distortion




subroutine get_local_asprat(du,dv,srat)
!
!
!  This subroutine computes the ratio of the square of the singular
!  values of the 3x2 jacobian matrix whose columsn are given by 
!  du and dv
!
!  J = [du dv]
!  J^{T}J = [|du|^2     <du,dv> ]
!           [<du,dv>      |dv|^2]
!
!  The singular values of J^{T}J are s1^2, s2^2 and are given by
!    (|du|^2 + |dv|^2 +- sqrt((|du|^2-|dv|^2)^2 + 4 <du,dv>^2))/2
!
!
!  input:
!    du - real *8 (3)
!       first column of jacobian matrix
!    dv - real *8 (3)
!       second column of jacobian matrix
!
!   output:
!      ratio: (s1/s2)^2
!


  implicit none
  real *8 du(3),dv(3),srat,a11,a12,a22,t,s1,s2

  a11 = du(1)**2 + du(2)**2 + du(3)**2
  a22 = dv(1)**2 + dv(2)**2 + dv(3)**2
  a12 = du(1)*dv(1) + du(2)*dv(2) + du(3)*dv(3)

  t = sqrt((a11-a22)**2+4*(a12**2))
  srat = (a11+a22+t)/(a11+a22-t)

end subroutine get_local_asprat








subroutine get_nximat(npatches,ixyzs,ixyzso,nximat)
  implicit none
  integer npatches,ixyzs(npatches+1),ixyzso(npatches+1),nximat,i

  nximat = 0
!$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:nximat) 
  do i=1,npatches
    nximat = nximat+(ixyzso(i+1)-ixyzso(i))*(ixyzs(i+1)-ixyzs(i))
  enddo
!$OMP END PARALLEL DO  

end subroutine get_nximat




subroutine get_ximats(npatches,iptype,norders,ixyzs,novers,ixyzso, &
  nximat,ximats,iximat)
  
  implicit none
  integer npatches,iptype(npatches),norders(npatches),ixyzs(npatches+1)
  integer novers(npatches),ixyzso(npatches+1),iximat(npatches),nximat
  integer, allocatable :: nximats(:)
  integer i,j,k,l,npols,kpols
  
  real *8 ximats(nximat)

  allocate(nximats(npatches))

!$OMP PARALLEL DO DEFAULT(SHARED)
  do i=1,npatches
    nximats(i) = (ixyzso(i+1)-ixyzso(i))*(ixyzs(i+1)-ixyzs(i))
  enddo
!$OMP END PARALLEL DO
 
 iximat(1) = 1
 nximats(1) = nximats(1)+1
 call cumsum(npatches-1,nximats,iximat(2))
 


!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(kpols,npols)

  do i=1,npatches
    if(iptype(i).eq.1) then
       kpols = ixyzs(i+1)-ixyzs(i)
       npols = ixyzso(i+1)-ixyzso(i)
       call koorn_oversamp_mat(norders(i),kpols,novers(i),npols, &
         ximats(iximat(i)))
    endif
  enddo
!$OMP END PARALLEL DO


end subroutine get_ximats
!
!
!
!
!
subroutine get_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffform)
!
!------------------------
!  This subroutine computes the first fundamental form at
!  the discretization nodes on the surface.
!
!  The first fundamental form is
!  
!  .. math::
!    
!    \begin{bmatrix} x_{u} \cdot x_{u} & x_{u} \cdot x_{v} \\
!    x_{u} \cdot x_{v} & x_{v} \cdot x_{v} \end{bmatrix}
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!  Output arguments:
!
!    - ffform: double precision(2,2,npts)
!        first fundamental form at the discretization nodes
!--------------------------
!
  
  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
  real *8, intent(out) :: ffform(2,2,npts)
  integer i
  real *8 xuu,xvv,xuv

  do i=1,npts
    xuu = srcvals(4,i)**2 + srcvals(5,i)**2 + srcvals(6,i)**2
    xvv = srcvals(7,i)**2 + srcvals(8,i)**2 + srcvals(9,i)**2
    xuv = srcvals(4,i)*srcvals(7,i) + srcvals(5,i)*srcvals(8,i)+ &
       srcvals(6,i)*srcvals(9,i)
    ffform(1,1,i) = xuu
    ffform(2,1,i) = xuv
    ffform(1,2,i) = xuv
    ffform(2,2,i) = xvv
  enddo

  return
end subroutine get_first_fundamental_form
!
!
!
!
!

subroutine get_inv_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffforminv)
!------------------------
!  This subroutine computes the inverse of the first fundamental form at
!  the discretization nodes on the surface.
!
!  The first fundamental form is
!  
!  .. math::
!    
!    \begin{bmatrix} x_{u} \cdot x_{u} & x_{u} \cdot x_{v} \\
!    x_{u} \cdot x_{v} & x_{v} \cdot x_{v} \end{bmatrix}
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!  Output arguments:
!
!    - ffforminv: double precision(2,2,npts)
!        inverse of first fundamental form at the discretization nodes
!--------------------------
  
  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches),npts
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
  real *8, intent(out) :: ffforminv(2,2,npts)
  integer i
  real *8 xuu,xvv,xuv,d

  do i=1,npts
    xuu = srcvals(4,i)**2 + srcvals(5,i)**2 + srcvals(6,i)**2
    xvv = srcvals(7,i)**2 + srcvals(8,i)**2 + srcvals(9,i)**2
    xuv = srcvals(4,i)*srcvals(7,i) + srcvals(5,i)*srcvals(8,i)+ &
       srcvals(6,i)*srcvals(9,i)
    d = xuu*xvv - xuv*xuv
    ffforminv(1,1,i) = xvv/d
    ffforminv(2,1,i) = -xuv/d
    ffforminv(1,2,i) = -xuv/d
    ffforminv(2,2,i) = xuu/d
  enddo

  return
end subroutine get_inv_first_fundamental_form
!
!
!
!
!
subroutine get_surf_grad(nd,npatches,norders,ixyzs,iptype,npts, &
  srccoefs,srcvals,f,df)
!
!-----------------------------
!  Compute the surface gradient of scalar function
!
!  Input arguments:
!
!    - nd: integer
!        number of functions
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!    - f: double precision (nd,npts)
!        function values at the discretization points
!
!  Output arguments:
!
!    - df: double precision(nd,2,npts)
!        surface gradient of the input functions
!-----------------------------
!
!

  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts,nd
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),f(nd,npts)
  real *8, intent(out) :: df(nd,2,npts)
  real *8, allocatable :: ffforminv(:,:,:)

  allocate(ffforminv(2,2,npts))


  call get_inv_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffforminv)

  call get_surf_grad_fast(nd,npatches,norders,ixyzs,iptype,npts, &
    ffforminv,f,df)

  return
end subroutine get_surf_grad
!
!
!
!
!
!
subroutine get_surf_grad_fast(nd,npatches,norders,ixyzs,iptype,npts, &
  ffforminv,f,df)
!
!-----------------------------
!  Faster version of computing the surface gradient of scalar function,
!  with precomputed inverse of first fundamental form
!
!  Input arguments:
!
!    - nd: integer
!        number of functions
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!    - npts: integer
!        total number of points on the surface
!    - ffforminv: double precision(2,2,npts)
!        inverse of first fundamental form at the discretization nodes
!    - f: double precision (nd,npts)
!        function values at the discretization points
!
!  Output arguments:
!
!    - df: double precision(nd,2,npts)
!        surface gradient of the input functions
!-----------------------------
!
!
  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts,nd
  real *8, intent(in) :: ffforminv(2,2,npts),f(nd,npts)
  real *8, intent(out) :: df(nd,2,npts)
  real *8, allocatable :: fcoefs(:,:),dfuv(:,:,:)
  integer i,idim

  allocate(dfuv(nd,2,npts))
  call get_surf_uv_grad(nd,npatches,norders,ixyzs,iptype,npts,f,dfuv)


  do i=1,npts
    do idim=1,nd
      df(idim,1,i) = dfuv(idim,1,i)*ffforminv(1,1,i) + &
         dfuv(idim,2,i)*ffforminv(1,2,i)
      df(idim,2,i) = dfuv(idim,1,i)*ffforminv(2,1,i) + &
         dfuv(idim,2,i)*ffforminv(2,2,i)
    enddo
  enddo


  return
end subroutine get_surf_grad_fast
!
!
!
!
!
subroutine get_surf_uv_grad(nd,npatches,norders,ixyzs,iptype,npts,f, &
   dfuv)
!
!-----------------------------
!  Compute the uv gradient of scalar function
!
!  Input arguments:
!
!    - nd: integer
!        number of functions
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!    - npts: integer
!        total number of points on the surface
!    - f: double precision (nd,npts)
!        function values at the discretization points
!
!  Output arguments:
!
!    - dfuv: double precision(nd,2,npts)
!        uv surface gradient of the input functions
!-----------------------------
!
  implicit none
  integer, intent(in) :: nd,npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1)
  integer, intent(in) :: iptype(npatches),npts
  real *8, intent(in) :: f(nd,npts)
  real *8, intent(out) :: dfuv(nd,2,npts)

  integer i,istart,npols



  do i=1,npatches
    istart = ixyzs(i)
    npols = ixyzs(i+1)-ixyzs(i)
    call get_surf_uv_grad_guru(nd,norders(i),npols,iptype(i), &
      f(1,istart),dfuv(1,1,istart))
  enddo


  return
end subroutine get_surf_uv_grad


subroutine get_surf_uv_grad_guru(nd,norder,npols,iptype,f,dfuv)
!
!-----------------------------
!  Compute the uv gradient of scalar function on a single triangular
!  patch
!
!  Input arguments:
!
!    - nd: integer
!        number of functions
!    - norder: integer
!        order of discretization
!    - npols: integer
!        total number of points on the patch
!    - f: double precision (nd,npols)
!        function values at the discretization points
!
!  Output arguments:
!
!    - dfuv: double precision(nd,2,npols)
!        uv surface gradient of the input functions
!-----------------------------
!
  implicit none
  integer, intent(in) :: nd,norder,npols,iptype
  real *8, intent(in) :: f(nd,npols)
  real *8, intent(out) :: dfuv(nd,2,npols)
  real *8 fcoefs(nd,npols),pols(npols),ders(2,npols)
  real *8 uv(2,npols)
  integer i,idim,j

  call vals_to_coefs_guru(nd,norder,npols,iptype,f,fcoefs)
  call get_disc_nodes(norder,npols,iptype,uv)

  do i=1,npols
    do idim=1,nd
      dfuv(idim,1,i) = 0
      dfuv(idim,2,i) = 0
    enddo
  enddo


  do i=1,npols
    call get_basis_ders(uv(1,i),norder,npols,iptype,pols,ders)
    do j=1,npols
      do idim=1,nd
        dfuv(idim,1,i) = dfuv(idim,1,i) + ders(1,j)*fcoefs(idim,j)
        dfuv(idim,2,i) = dfuv(idim,2,i) + ders(2,j)*fcoefs(idim,j)
      enddo
    enddo
  enddo

  return
end subroutine get_surf_uv_grad_guru
!
!
!
!
!
subroutine get_surf_div(nd,npatches,norders,ixyzs,iptype,npts, &
  srccoefs,srcvals,f,df)
!
!-----------------------------
!  Compute the surface divergence of a vector function
!
!  Input arguments:
!
!    - nd: integer
!        number of vector fields
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!    - f: double precision (nd,2,npts)
!        vector field values at the discretization points
!
!  Output arguments:
!
!    - df: double precision(nd,npts)
!        surface divergence of the input vector fields
!-----------------------------
!
  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts,nd
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
  real *8, intent(in) :: f(nd,2,npts)
  real *8, intent(out) :: df(nd,npts)
  real *8, allocatable :: ffform(:,:,:)

  allocate(ffform(2,2,npts))


  call get_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffform)

  call get_surf_div_fast(nd,npatches,norders,ixyzs,iptype,npts, &
    ffform,f,df)

  return
end subroutine get_surf_div
!
!
!
!
!


subroutine get_surf_div_fast(nd,npatches,norders,ixyzs,iptype,npts, &
  ffform,f,df)
!
!-----------------------------
!  Faster version of computing the surface divergence of a vector function
!  with precmoputed first fundamental form
!
!  Input arguments:
!
!    - nd: integer
!        number of vector fields
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!    - npts: integer
!        total number of points on the surface
!    - f: double precision (nd,2,npts)
!        vector field values at the discretization points
!
!  Output arguments:
!
!    - df: double precision(nd,npts)
!        surface divergence of the input vector fields
!-----------------------------
!
  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts,nd
  real *8, intent(in) :: ffform(2,2,npts),f(nd,2,npts)
  real *8, intent(out) :: df(nd,npts)
  integer i,istart,npols

  do i=1,npatches
    istart = ixyzs(i)
    npols = ixyzs(i+1)-ixyzs(i)
    call get_surf_div_guru(nd,norders(i),npols,iptype(i), &
      ffform(1,1,istart),f(1,1,istart),df(1,istart))
  enddo


  return

  return
end subroutine get_surf_div_fast
!
!
!
!
!
subroutine get_surf_div_guru(nd,norder,npols,iptype,ff,f,df)
!
!-----------------------------
!  Compute the surface divergence of a vector field on a single triangular
!  patch
!
!  Input arguments:
!
!    - nd: integer
!        number of functions
!    - norder: integer
!        order of discretization
!    - npols: integer
!        total number of points on the patch
!    - ff: double precision(2,2,npols)
!        first fundamental form at discretization nodes
!    - f: double precision (nd,2,npols)
!        vector field values at the discretization points
!
!  Output arguments:
!
!    - df: double precision(nd,npols)
!        surface divergecen of the input vector fields 
!-----------------------------
!
  implicit none
  integer, intent(in) :: nd,norder,npols,iptype
  real *8, intent(in) :: f(nd,2,npols),ff(2,2,npols)
  real *8, intent(out) :: df(nd,npols)
  real *8 fcoefs(nd,2,npols),pols(npols),ders(2,npols)
  real *8 ftmp(nd,2,npols)
  real *8 uv(2,npols),gg(npols)
  integer i,idim,j

  do i=1,npols
    gg(i) = sqrt(ff(1,1,i)*ff(2,2,i) - ff(1,2,i)*ff(2,1,i))
    do idim=1,nd
      ftmp(idim,1,i) = f(idim,1,i)*gg(i)
      ftmp(idim,2,i) = f(idim,2,i)*gg(i)
    enddo
  enddo

  call vals_to_coefs_guru(2*nd,norder,npols,iptype,ftmp,fcoefs)
  call get_disc_nodes(norder,npols,iptype,uv)

  do i=1,npols
    do idim=1,nd
      df(idim,i) = 0
    enddo
  enddo


  do i=1,npols
    call get_basis_ders(uv(1,i),norder,npols,iptype,pols,ders)
    do j=1,npols
      do idim=1,nd
        df(idim,i) = df(idim,i) + ders(1,j)*fcoefs(idim,1,j) + &
                                  ders(2,j)*fcoefs(idim,2,j)
      enddo
    enddo

    do idim=1,nd
      df(idim,i) = df(idim,i)/gg(i)
    enddo
  enddo

  return
end subroutine get_surf_div_guru
!
!
!
!

subroutine col_ind_to_patch_node_ind(npatches,ixyzs,ncol,col_ind, &
  patch_ind,node_ind)

  implicit none
  integer, intent(in) :: npatches,ixyzs(npatches+1),ncol
  integer, intent(in) :: col_ind(ncol)
  integer, intent(out) :: patch_ind(ncol),node_ind(ncol)

  integer i,j,j0,npols

  npols = ixyzs(2)-ixyzs(1)
  do i=1,ncol
    node_ind(i) = mod(col_ind(i),npols)
    patch_ind(i) = col_ind(i)/npols+1
    if(node_ind(i).eq.0) then
      node_ind(i) = npols
      patch_ind(i) = patch_ind(i) - 1
    endif
  enddo

!  call sorti(ncol,col_ind,iper)
!
!  do i=1,ncol
!    col_ind_sort(i) = col_ind(iper(i))
!    patch_ind(i) = 0
!    node_ind(i) = 0
!  enddo
!
!
!  j0 = 1
!  do i=1,ncol
!    do j=j0,npatches
!      if(ixyzs(j+1).gt.col_ind_sort(i)) then
!        patch_ind(iper(i)) = j
!        node_ind(iper(i)) = col_ind(iper(i))-ixyzs(j)+1
!        j0 = j
!        goto 1111
!      endif
!    enddo
! 1111 continue  
!  enddo
  

end subroutine col_ind_to_patch_node_ind




subroutine plot_surface_info_all(dlam,npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,fname,title)
!
!-----------------------------
!  This subroutine creates a vtk file which contains the following scalar fields:
!    * points per wavelength
!    * patch distortion
!    * estimated error in surface representation if order>=2, else error is
!      set to 0
!
!  Input arguments:
!    - dlam: real *8
!        wavelength
!    - nd: integer
!        number of functions
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!    - fname: character (len=*)
!        file name to which info needs to be written
!
   implicit none 
   real *8, intent(in) :: dlam
   integer, intent(in) :: npatches,norders(npatches),ixyzs(npatches+1)
   integer, intent(in) :: iptype(npatches),npts
   real *8, intent(in) :: srcvals(12,npts),srccoefs(9,npts)
   character (len=*) :: fname,title

   real *8, allocatable :: dppw(:),pdis(:),errp(:)
   real *8, allocatable :: sigma(:,:),cms(:,:),rads(:),qwts(:)
   character (len=10), dimension(3) :: scalar_names

   integer nsc,nd,ip,j,i
   real *8 errm
   

   nsc = 10

   scalar_names(1) = 'ppw'
   scalar_names(2) = 'distortion'
   scalar_names(3) = 'surf_err'

   allocate(dppw(npatches),pdis(npatches),errp(npatches))
   allocate(cms(3,npatches),rads(npatches))
   call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, &
    srccoefs,cms,rads)

   do i=1,npatches
     dppw(i) = (norders(i)+0.0d0)/(rads(i)/dlam)
   enddo

   allocate(qwts(npts))
   call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,qwts)


   call get_patch_distortion(npatches,norders,ixyzs,iptype,npts,&
    srccoefs,srcvals,qwts,pdis)

   call surf_fun_error(9,npatches,norders,ixyzs,iptype, &
     npts,srcvals(1:9,1:npts),qwts,errp,errm)

   nd = 3
   allocate(sigma(3,npts))
   do ip=1,npatches
     do j=ixyzs(ip),ixyzs(ip+1)-1
       sigma(1,j) = dppw(ip)
       sigma(2,j) = pdis(ip)
       sigma(3,j) = errp(ip)
     enddo
   enddo

   call surf_vtk_plot_scalar_many(nd,npatches,norders,ixyzs,iptype, &
     npts,srccoefs,srcvals,sigma,nsc,scalar_names,fname,title)



end subroutine plot_surface_info_all
!
!
!
!
!
!


subroutine getpatchinfo(npatches, patchpnt, par1, par2, par3, par4, &
    npols, uvs, umatr, srcvals, srccoefs)
  implicit real *8 (a-h,o-z)
  real *8 :: uvs(2,npols), srcvals(12,*)
  real *8 :: umatr(npols,npols),srccoefs(9,*)
  external patchpnt

  real *8 :: xyz(3), dxyzduv(3,10), xyznorm(3)
  real *8 :: xyztang1(3), xyztang2(3)

  !
  !       This subroutine return all points, normals and tangents from
  !       geometry descriptor
  !
  !       Input parameters:
  !
  !         npatches: integer: the number of patches
  !         patchpnt: external: subroutine evaluating points along
  !               the surface, given patch by patch, must be of the form
  !                     patchpnt(ipatch,u,v,xyz,dxyzduv,par1,par2,par3,par4)
  !         par1,par2,par3,par4: extra parameters
  !         npols: integer: the total number of polynomials for each patch
  !         uvs: real *8(2,npols): local u-discretization points for each patch
  !         umatr: real *8(npols,npols): values to coeffs mapatchx on standard patch 
  !
  !       Output parameters:
  !
  !         srcvals: real *8(12,npts): geometry info with first derivatives
  !               srcvals(1:3,:) - xyz
  !               srcvals(4:6,:) - dxyz/du
  !               srcvals(7:9,:) - dxyz/dv
  !               srcvals(10:12,:) - xyznorms
  !
  !         srccoefs: real *8 (9,npts): geometry info as koornwinder expansion
  !                     coefficients
  !                    
  !         npts: integer: the total number of points in discretization
  !


  do ipatch=1,npatches
    do i=1,npols

      u=uvs(1,i)
      v=uvs(2,i)

      ipt = (ipatch-1)*npols + i 

      call patchpnt(ipatch,u,v,srcvals(1,ipt),srcvals(4,ipt),par1, &
             par2,par3,par4)

      call cross_prod3d(srcvals(4,ipt),srcvals(7,ipt),srcvals(10,ipt))

      ds = sqrt(srcvals(10,ipt)**2 + srcvals(11,ipt)**2 + &
              srcvals(12,ipt)**2)
      srcvals(10,ipt) = srcvals(10,ipt)/ds
      srcvals(11,ipt) = srcvals(11,ipt)/ds
      srcvals(12,ipt) = srcvals(12,ipt)/ds

    end do

    do i=1,npols
      ipt = (ipatch-1)*npols + i
      do j=1,9
        srccoefs(j,ipt) = 0
        do l=1,npols
          lpt = (ipatch-1)*npols + l
          srccoefs(j,ipt) = srccoefs(j,ipt) + umatr(i,l)*srcvals(j,lpt)
        end do
      end do
    end do
  end do

  npts = npatches*npols

  return
end subroutine getpatchinfo


subroutine get_surf_interp_mat_targ_mem(npatches,ixyzs,ntarg, &
  ipatchtarg,lmem)
!
!  This subroutine estimates the memory required for storing an
!  interpolation matrix for a given collection of targets
!  on surface whose patch id have already been identified
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - ntarg: integer
!        number of targets
!    - ipatchtarg: integer(ntarg)
!        patch id of target i
!
!  Output arguments:
!    - lmem: integer
!        memory required for storing the interpolation matrix

  implicit none
  integer, intent(in) :: npatches,ixyzs(npatches+1),ntarg
  integer, intent(in) :: ipatchtarg(ntarg)
  integer, intent(out) :: lmem

!  Temporary variables

  integer i,ipatch,itarg

  lmem = 0
!$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:lmem) private(i) 
  do itarg=1,ntarg
    i = ipatchtarg(itarg)
    lmem = lmem + ixyzs(i+1) - ixyzs(i) 
  enddo
!$OMP END PARALLEL DO  

end subroutine get_surf_interp_mat_targ_mem
!
!
!
!
!
subroutine get_surf_interp_mat_targ(npatches,norders,ixyzs,iptype,npts,ntarg, &
  ipatchtarg,uvs_targ,lmem,xmattarg,ixmattarg)
!
!  This subroutine constructs the
!  interpolation matrix for a given collection of targets
!  on surface whose patch id have already been identified
!
!  This subroutine currently implements a slow version
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        discretization order of patches
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!    - npts: integer
!        number of points on surface
!    - ntarg: integer
!        number of targets
!    - ipatchtarg: integer(ntarg)
!        patch id of target i
!    - uvs_targ: real *8 (2,ntarg)
!        local uv coordinates of targets on surface
!    - lmem: integer
!        memory required to store the interpolation matrix
!
!  Output arguments:
!    - xmattarg: real *8 (lmem)
!        the interpolation matrix
!    - ixmattarg: integer(ntarg+1)
!        location where interpolation matrix for target i begins
  implicit none
  integer, intent(in) :: npatches,npts
  integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
  integer, intent(in) :: iptype(npatches)
  integer, intent(in) :: ntarg
  integer, intent(in) :: ipatchtarg(ntarg)
  real *8, intent(in) :: uvs_targ(2,ntarg)
  integer, intent(in) :: lmem
  real *8, intent(out) :: xmattarg(lmem)
  integer, intent(out) :: ixmattarg(ntarg+1)

!  Temporary variables

  integer i,ipatch,itarg,lmem0
  integer, allocatable :: nxmattarg(:)

  real *8, allocatable :: umat(:,:,:),pols(:),uvs(:,:)
  integer, allocatable :: iuni(:),iuniind(:),nordertarg(:)
  integer incx,incy,iind,norder,npmax,npols,nuni
  real *8 alpha,beta

  allocate(nxmattarg(ntarg))

!$OMP PARALLEL DO DEFAULT(SHARED) private(i) 
  do itarg=1,ntarg
    i = ipatchtarg(itarg)
    nxmattarg(itarg) = ixyzs(i+1) - ixyzs(i) 
  enddo
!$OMP END PARALLEL DO  


 ixmattarg(1) = 1
 nxmattarg(1) = nxmattarg(1)+1
 call cumsum(ntarg,nxmattarg,ixmattarg(2))

 npmax = maxval(nxmattarg)
 allocate(pols(npmax))

 allocate(nordertarg(ntarg))

!$OMP PARALLEL DO DEFAULT(SHARED)
 do i=1,ntarg
   nordertarg(i) = norders(ipatchtarg(i))
 enddo
!$OMP END PARALLEL DO 


 allocate(iuni(ntarg),iuniind(ntarg))
 call get_iuni1_omp(ntarg,nordertarg,nuni,iuni,iuniind) 
 allocate(umat(npmax,npmax,nuni),uvs(2,npmax))

 do i=1,nuni
   norder = iuni(i)
   if(iptype(iuni(i)).eq.1) then
     npols = (norder+1)*(norder+2)/2
     call get_vioreanu_nodes(norder,npols,uvs)
     call koorn_vals2coefs(norder,npols,uvs,umat(1:npols,1:npols,i))
   endif
 enddo

 alpha = 1.0d0
 beta = 0.0d0
 incx = 1
 incy = 1

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(itarg,ipatch,iind,pols,npols,norder)
  do itarg=1,ntarg
    ipatch = ipatchtarg(itarg)
    iind = iuniind(itarg)
    norder = norders(ipatch)
    npols = (norder+1)*(norder+2)/2 
    call koorn_pols(uvs_targ(1,itarg),norders(ipatch),npols,pols)

    call dgemv_guru('t',npols,npols,alpha,umat(1,1,iind),npmax,pols,incx, &
      beta,xmattarg(ixmattarg(itarg)),incy)

  enddo
!$OMP END PARALLEL DO 

end subroutine get_surf_interp_mat_targ
!
!
!
!
!
subroutine get_second_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,sfform)
!
!------------------------
!  This subroutine computes the second fundamental form at
!  the discretization nodes on the surface.
!
!  The second fundamental form is
!  
!  .. math::
!    
!    \begin{bmatrix} x_{uu} \cdot n & x_{uv} \cdot n \\
!    x_{uv} \cdot n & x_{vv} \cdot n \end{bmatrix}
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!  Output arguments:
!
!    - sfform: double precision(2,2,npts)
!        second fundamental form at the discretization nodes
!--------------------------
!
  
  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
  real *8, intent(out) :: sfform(2,2,npts)
  integer i
  real *8 L,M,N
  real *8, allocatable :: dxuv(:,:)
  real *8, allocatable :: dxuv2(:,:,:)


  allocate(dxuv(6,npts))
! Calculating x_{uu}, x_{uv}, x_{uv} stored in dxuv2 
 
  do i=1,npts
    dxuv(1,i) = srcvals(4,i)  
    dxuv(2,i) = srcvals(5,i)  
    dxuv(3,i) = srcvals(6,i)  
    dxuv(4,i) = srcvals(7,i)  
    dxuv(5,i) = srcvals(8,i)  
    dxuv(6,i) = srcvals(9,i)   
  enddo

  allocate(dxuv2(6,2,npts))

  call get_surf_uv_grad(6,npatches,norders,ixyzs,iptype,npts,dxuv,dxuv2)



  do i=1,npts
    L = dxuv2(1,1,i)*srcvals(10,i) + dxuv2(2,1,i)*srcvals(11,i) + dxuv2(3,1,i)*srcvals(12,i) ! Calculation of L, M, N. L = x_uu \cdot n 
    M = dxuv2(1,2,i)*srcvals(10,i) + dxuv2(2,2,i)*srcvals(11,i) + dxuv2(3,2,i)*srcvals(12,i)  ! M = x_uv \cdot n
    N = dxuv2(4,2,i)*srcvals(10,i) + dxuv2(5,2,i)*srcvals(11,i) + dxuv2(6,2,i)*srcvals(12,i)  ! N = x_vv \cdot n
    sfform(1,1,i) = L
    sfform(2,1,i) = M
    sfform(1,2,i) = M
    sfform(2,2,i) = N
  enddo

  return
end subroutine get_second_fundamental_form
!
!
!
subroutine get_mean_curvature(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,mean_curv)
!
!------------------------
!  This subroutine computes the mean curvature at
!  the discretization nodes on the surface.
!
!  The mean curvature is
!  
!  .. math::
!    
!    0.5*Trace(II \cdot I^{-1}) \\
!    
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!       * iptype = 1, triangular patch with RV nodes
!       * iptype = 11, quad patch with GL nodes
!       * iptype = 12, quad patch with Cheb nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!  Output arguments:
!
!    - mean_curv: double precision(npts)
!        mean curvature at the discretization nodes
!--------------------------
!
  
  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
  real *8, intent(out) :: mean_curv(npts)
  integer i
  real *8, allocatable :: ffform(:,:,:),ffforminv(:,:,:)
  real *8, allocatable :: sfform(:,:,:)

  allocate(ffform(2,2,npts))

  call get_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffform)



  allocate(ffforminv(2,2,npts))

  call get_inv_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffforminv)

  allocate(sfform(2,2,npts))

  call get_second_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,sfform)
  
!  print *,"Point on surface:", srcvals(1,3),srcvals(2,3), srcvals(3,3) 
!  print *,"First fundamental form=", ffform(:,:, 3) 
!  print *,"Inverse first fundamental form=", ffforminv(:,:, 3)
!  print *,"Second fundamental form=", sfform(:,:, 3)
 


! Calculating mean curvature 
 
  do i=1,npts
    mean_curv(i) = -0.5*(sfform(1,1,i)*ffforminv(1,1,i) + &
                     sfform(1,2,i)*ffforminv(2,1,i) + &
                     sfform(2,1,i)*ffforminv(1,2,i) + &
                     sfform(2,2,i)*ffforminv(2,2,i))
  enddo
!  print *,"Mean=", mean_curv(3)
 
  return
end subroutine get_mean_curvature

!
!
!
