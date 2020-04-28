!  This file contains the following user callable routines
!
!
!      get_far_order - determine oversampling order for the far-field
!          based on bounding radii for near-field, the accuracy and the kernel
!
!      get_far_order_tri - called by get_far_order when patch is a triangle
!
!
subroutine get_far_order(eps,npatches,norders,ixyzs,iptype,cms,rads,&
  npts,srccoefs, ndtarg, ntarg,targvals,ikerorder,zk,nnz,row_ptr, &
  col_ind,rfac,nfars,ixyzso)

!
!
!  For a given surface, this subroutine computes the 
!  quadrature order required to compute the far-field accurately
!
!  The method proceeds by computing for each patch,
!  the potential e^{ikr}/r (if ikerorder=-1)
!  \nabla e^{ikr}/r \cdot n if ikerorder=0 and
!  \Hessian e^{ikr}/r \cdot n \cdot n if ikerorder=1
!  at the 10 furthest
!  targets, and 15 targets randomly sprinkled on the boundary of
!  the defintion of the near-field using a RV nodes of order 1,2,3,4,6,8,10,12
!  and stops whenever the potential at all the targets has achieved the
!  desired precision (in absolute value)
!
!  input
!    eps - real *8
!       tolerance requested
!
!    npatches - integer
!       number of patches
!
!    norders - integer
!       order of discretization for each patch
!
!    ixyzs - integer(npatches+1)
!       ixyzs(i) indicates the location in srccoefs array
!       where information for patch i begins
!
!    npts - integer
!       total number of discretization points on the surface
!
!    srccoefs - real *8 (9,npts)
!       orthogonal polynomial expansions on each patch
!       of x,y,z,dxyz/du,dxyz/dv
!
!    ndtarg - integer
!       leading dimension of target array (the first three
!       dimensions must be target locations)
!    
!    ntarg - integer
!       number of targets
!
!    targvals - real *8 (ndtarg,ntarg)
!       target info in the near field, note the first three 
!       inputs per point must be the xyz coordinates
!
!    ikerorder - integer
!       type of kernel
!         ikerorder = -1 -> single layer type operator
!         ikerorder = 0 -> Double layer type operator
!         ikerorder = 1 -> derivative of double layer
!
!    zk - complex *16
!         Wave number for the problem
!    
!    nnz - integer
!       number of non-zero interactions in the near-field
!
!    row_ptr - integer (ntarg+1)
!       row_ptr(i) indicates the location in col_ind array
!       where list of source patches relevant to target i
!       start
!
!    col_ind - integer(nnz)
!       list of source patch indices in the near field
!       of targets. col_ptr(row_ind(i):row_ind(i+1)-1))
!       is the collection of source patches relevant 
!       for target i
!
!    rfac - real *8
!       factor used for determining the near-field 
!       for each patch
!    
!    output
!      nfars - integer(npatches)
!        quadrature order for far-field for each patch
!
!      ixyzso - integer(npatches+1)
!         location in oversampled source array where
!         information of patch i starts
!
!
  implicit real *8 (a-h,o-z)
  integer npatches,norders(npatches),ixyzs(npatches),iptype(npatches)
  integer npts,ikerorder
  real *8 cms(3,npatches),rads(npatches),srccoefs(9,npts)
  real *8 dpars
  complex *16 zk
  integer ipars

  integer ndtarg,ntarg
  real *8 targvals(ndtarg,ntarg)
  integer nnz
  integer row_ptr(ntarg+1),col_ind(nnz)
  real *8 rfac
  integer nfars(npatches),ixyzso(npatches+1)

  integer, allocatable :: col_ptr(:),row_ind(:),iper(:)

  real *8, allocatable :: targtmp(:,:)
  real *8, allocatable :: uvs(:,:),wts(:),umat(:,:),vmat(:,:)
  character *1 transa,transb

  allocate(col_ptr(npatches+1),row_ind(nnz),iper(nnz))

  call rsc_to_csc(npatches,ntarg,nnz,row_ptr,col_ind,col_ptr,&
    row_ind,iper)

  ixyzso(1) = 1
  do i=1,npatches
    npols = ixyzs(i+1)-ixyzs(i)
    istart = ixyzs(i)
    if(iptype(i).eq.1) then
       norder = norders(i)


!
!      gather the targets
!
      ntarg0 = col_ptr(i+1)-col_ptr(i)
      allocate(targtmp(ndtarg,ntarg0))

      ii = 0
      do j=col_ptr(i),col_ptr(i+1)-1
        jpt = row_ind(j)
        ii = ii +1
        do l=1,ndtarg
          targtmp(l,ii) = targvals(l,jpt)
        enddo
      enddo

      rad0 = rfac*rads(i)



      
      call get_far_order_tri(eps,norder,npols,cms(1,i),rad0, &
        srccoefs(1,istart),ikerorder,zk, &
        ndtarg,ntarg0,targtmp,nfars(i),npolsf)
      ixyzso(i+1) = ixyzso(i)+npolsf
      deallocate(targtmp)

    endif
  enddo
end subroutine get_far_order





subroutine get_far_order_tri(eps,norder,npols,cm,rad,srccoefs, &
  ikerorder,zk,ndtarg,ntarg,targvals,nfar,npolsf)

!
!  NOTE: This routine is not currently optimized for performance
!
!
!  For a given triangular patch, this subroutine computes the 
!  quadrature order required to compute the far-field accurately
!
!  The method proceeds by computing the potential e^{ikr}/r (if ikerorder=-1)
!  \nabla e^{ikr}/r \cdot n if ikerorder=0 and
!  \Hessian e^{ikr}/r \cdot n \cdot n if ikerorder=1
!  at the 10 furthest
!  targets, and 15 targets randomly sprinkled on the boundary of
!  the defintion of the near-field using a RV nodes of order 1,2,3,4,6,8,10,12
!  and stops whenever the potential at all the targets has achieved the
!  desired precision (in absolute value)
!
!  input
!    eps - real *8
!       tolerance requested
!
!    norder - integer
!       order of patch
!
!    npols - integer
!       number of discretization nodes on patch
!       npols = (norder+1)*(norder+2)/2
!
!    cm - real *8 (3)
!       centroid of patch
!
!    rad - real *8 
!       radius of definition of near-field
!       far-field assumed to be defined via |x-cm|>rad
!
!    srccoefs - real *8 (9,npols)
!       koornwinder expansion coeffs
!
!    ikerorder - integer
!       type of kernel
!         ikerorder = -1 -> single layer type operator
!         ikerorder = 0 -> Double layer type operator
!         ikerorder = 1 -> derivative of double layer
!
!    zk - complex *16
!         Wave number for the problem
!
!    ndtarg - integer
!       leading dimension of target array (the first three
!       dimensions must be target locations)
!    
!    ntarg - integer
!       number of targets
!
!    targvals - real *8 (ndtarg,ntarg)
!       target info in the near field, note the first three 
!       inputs per point must be the xyz coordinates
!    
!    output
!      nfar - integer
!        quadrature order for far-field (currently will only be 1,2,3,4,6,8,10,12)
!
!
  implicit real *8 (a-h,o-z)
  real *8 eps,cm(3),rad,srccoefs(3,npols)
  real *8 targvals(ndtarg,ntarg)
  integer norder,npols,ndtarg,ntarg,nfar
  complex *16 zk
  real *8, allocatable :: targtest(:,:),dd(:)
  integer, allocatable :: indd(:)
  real *8 phi,thet,pi,done
  integer i,ii

  real *8 dpars
  integer ipars

  integer ndd,ndz,ndi,ndt

  real *8, allocatable ::  srctmp(:,:)
  complex *16, allocatable :: sigmatmp(:),vcomp(:,:),vref(:,:)
  complex *16 val
  real *8, allocatable :: pmat(:,:)
  real *8, allocatable :: uvs(:,:),wts(:),qwts(:)
  real *8 tmp(3)
  integer nfars(8)
  integer nmax,iseed1,iseed2 

  character *1 transa,transb

  procedure (), pointer :: fker

  external h3d_slp,h3d_dlp,h3d_qlp

  if(ikerorder.eq.-1) fker => h3d_slp
  if(ikerorder.eq.0) fker => h3d_dlp
  if(ikerorder.eq.1) fker => h3d_qlp


  ndd = 0
  ndz = 1
  ndi = 0
  ndt = 3

  iseed1 = 1001
  iseed2 = 3042

  nfars(1) = 1
  nfars(2) = 2
  nfars(3) = 3
  nfars(4) = 4
  nfars(5) = 6
  nfars(6) = 8
  nfars(7) = 10
  nfars(8) = 12

  done = 1
  pi = atan(done)*4

!
!  First extract relevant sets of targets
!
   allocate(dd(ntarg),indd(ntarg))

!     max number of targets
   nmax = 25

   allocate(targtest(3,nmax))

   do i=1,ntarg
     dd(i) = (targvals(1,i)-cm(1))**2 + (targvals(2,i)-cm(2))**2 + &
       (targvals(3,i)-cm(3))**2
   enddo

   call sortr(ntarg,dd,indd)




   i0 = min(10,ntarg/2)

   do i=1,i0
     ii = indd(ntarg-i+1)

     do j=1,ndtarg
       targtest(j,i) = targvals(j,ii)
     enddo
   enddo

   do i=i0+1,nmax
     phi = hkrand(iseed1+i)*2*pi
     thet = hkrand(iseed2+i)*pi
     targtest(1,i) = cm(1) + rad*sin(thet)*cos(phi)
     targtest(2,i) = cm(2) + rad*sin(thet)*sin(phi)
     targtest(3,i) = cm(3) + rad*cos(thet)
   enddo


!
!  end of generating targets
!


  nomax = 13
  npmax = (nomax+1)*(nomax+2)/2
  allocate(srctmp(12,npmax),sigmatmp(npmax),pmat(npols,npmax))
  allocate(uvs(2,npmax),wts(npmax),qwts(npmax))

  allocate(vcomp(npols,nmax),vref(npols,nmax))


  iistart = 0
  do i=1,8
    if(nfars(i).lt.norder) iistart = i 
  enddo

  iistart = max(iistart,1)
  iistart = 1


  nfar = nfars(iistart) 
  npolsf = (nfar+1)*(nfar+2)/2
  call get_vioreanu_nodes_wts(nfar,npolsf,uvs,wts)


  do i=1,npolsf
    call koorn_pols(uvs(1,i),norder,npols,pmat(1,i)) 
  enddo


  transa = 'n'
  transb = 'n'
  alpha = 1
  beta = 0
  
  call dgemm(transa,transb,9,npolsf,npols,alpha,srccoefs,9,pmat,npols,&
     beta,srctmp,12)

!
!   compute the potential at the nmax targets
!
  do jpt=1,npolsf
    call cross_prod3d(srctmp(4,jpt),srctmp(7,jpt),tmp)
    rr = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
    srctmp(10,jpt) = tmp(1)/rr
    srctmp(11,jpt) = tmp(2)/rr
    srctmp(12,jpt) = tmp(3)/rr
    qwts(jpt) = wts(jpt)*rr
  enddo

  do itarg=1,nmax
    do l=1,npols
      vcomp(l,itarg) = 0
    enddo
    do j=1,npolsf
      call fker(srctmp(1,j),ndt,targtest(1,itarg),ndd,dpars,ndz,zk,ndi, &
        ipars,val)
      do l=1,npols
        vcomp(l,itarg) = vcomp(l,itarg)+val*pmat(l,j)*qwts(j)
      enddo
    enddo
  enddo




  do nnn=iistart+1,8
    nfar = nfars(nnn)

    npolsf = (nfar+1)*(nfar+2)/2
    call get_vioreanu_nodes_wts(nfar,npolsf,uvs,wts)

  
    do i=1,npolsf
      call koorn_pols(uvs(1,i),norder,npols,pmat(1,i)) 
    enddo


    transa = 'n'
    transb = 'n'
    alpha = 1
    beta = 0
  
    call dgemm(transa,transb,9,npolsf,npols,alpha,srccoefs,9,pmat, &
       npols,beta,srctmp,12)


!
!   compute the potential at the nmax targets
!   using order 4 nodes 
!
    do jpt=1,npolsf
      call cross_prod3d(srctmp(4,jpt),srctmp(7,jpt),tmp)
      rr = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
      srctmp(10,jpt) = tmp(1)/rr
      srctmp(11,jpt) = tmp(2)/rr
      srctmp(12,jpt) = tmp(3)/rr
      qwts(jpt) = wts(jpt)*rr
    enddo

    errmax=0 
    do itarg=1,nmax
      do l=1,npols
        vref(l,itarg) = 0
      enddo
      do j=1,npolsf
        call fker(srctmp(1,j),ndt,targtest(1,itarg),ndd,dpars,ndz,zk,ndi, &
        ipars,val)
        do l=1,npols
          vref(l,itarg) = vref(l,itarg)+val*pmat(l,j)*qwts(j)
        enddo
      enddo

      do l=1,npols
        err = abs(vref(l,itarg)-vcomp(l,itarg))
        if(err.gt.errmax) errmax = err
      enddo
    enddo



    if(errmax.lt.eps) goto 1111
    do itarg=1,nmax
      do l=1,npols
        vcomp(l,itarg) = vref(l,itarg)
      enddo
    enddo
  enddo

 1111 continue
  
  iistart = nnn-1
  iistart = max(iistart,1)
  nfar = nfars(iistart)
  npolsf = (nfar+1)*(nfar+2)/2 


end subroutine get_far_order_tri
!
!
!
