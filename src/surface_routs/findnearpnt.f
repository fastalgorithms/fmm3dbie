c
c   This file has two user callable routines.
c
c   * get_closest_points: Given a collection of targets, and a surface
c       for any target that is within a tubular neighborhood of the
c       surface, it identifies the closest point on surface
c
c   * findnear_surface_point: Given a collection of targets, a surface
c      and a closest discretization point on surface for each target,
c      it identifies the closest point on the patch
c
c
      subroutine get_closest_points(ndtarg, ntarg, targs,
     1     npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals,
     2     ipatch_id, uvs_src, maxiter, tol, rfac, sxyz, 
     3     ipatch_id_targ, uvsloc, dists, flags)
c
c  Given a collection of targets and a surface, this subroutine
c  identifies the closest point on the surface if the target lies
c  within a tubuluar neighborhood of the surface given by
c  {x | d(x,c_{j}) < rfac*r_{j}} where c_j is the centroid of patch
c  j and r_{j} is the radius of the bounding sphere.
c
c  Input arguments:
c    - ndtarg: integer *8
c        leading dimension of target array, should at least be 3.
c    - ntarg: integer *8
c        number of targets
c    - targs: real *8(ndtarg,ntarg)
c        The first three components must be the
c        xyz coordinates of the targets
c    - npatches: integer *8
c        number of patches
c    - norders: integer *8(npatches)
c        order of discretization on each patch
c    - ixyzs: integer *8(npatches+1)
c        ixyzs(i) denotes the starting location in srccoefs,
c        and srcvals array where information for patch i begins
c    - iptype: integer *8(npatches)
c        type of patch
c       * iptype = 1, triangular patch with RV nodes
c       * iptype = 11, quad patch with GL nodes
c       * iptype = 12, quad patch with Cheb nodes
c    - npts: integer *8
c        total number of discretization points on the boundary
c    - srccoefs: double precision (9,npts)
c        basis expansion coefficients of x, $\partial_{u} x$,
c        and $\partial_{v} x$.
c    - srcvals: double precision (12,npts)
c        x, $\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
c        discretization nodes
c    - ipatch_id: integer *8(npts)
c        patch id of the discretization points
c    - uvs_src: real *8(2,npts)
c        local uv coordinates of the discretization points
c    - maxiter: integer *8
c        maximum number of iterations to run newton for
c    - tol: real *8
c        tolerance to be used for the Newton algorithm
c    - rfac: real *8
c        factor defining the tubular neighborhood, see 
c        documentation above
c
c   Output arguments:
c    - sxyz: real *8 (3,ntarg)
c        xyz coordinates of the closest point on the patch 
c        corresponding to the target-patch input. Only
c        meaningful if flags(i) = 0
c    - ipatch_id_targ: integer *8(ntarg)
c        patch number of closest point on the patch.
c        Only meaningful if flags(i) = 0
c    - uvsloc: real *8(2,ntarg)
c        uv coordinates of the closest point on the patch.
c        Only meaningful if flags(i) = 0
c    - dists: real *8(ntarg)
c        distance of the closest point to the target
c        Only meaningful if flags(i) = 0
c    - flags: integer *8(ntarg)
c        flag corresponding to succesful execution of newton
c        flags(i) = 0 => Newton for target i converged successfully.
c        flags(i) = -1 => target not in tubular neighborhood
c        flags(i) = 1 => Newton failed execution
c-------------------------------------
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: ntarg, ndtarg
      real *8, intent(in) :: targs(ndtarg, ntarg)

      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srcvals(12,npts), srccoefs(9,npts)
      integer *8, intent(in) :: ipatch_id(npts)
      real *8, intent(in) :: uvs_src(2,npts)

      integer *8, intent(in) :: maxiter
      real *8, intent(in) :: tol
      
      real *8, intent(in) :: rfac

      real *8, intent(out) :: sxyz(3,ntarg)
      integer *8, intent(out) :: ipatch_id_targ(ntarg)
      real *8, intent(out) :: uvsloc(2,ntarg)
      real *8, intent(out) :: dists(ntarg)
      integer *8, intent(out) :: flags(ntarg)

      real *8, allocatable :: cms(:,:), rads(:)
      integer *8, allocatable :: row_ptr(:), col_ind(:)
      real *8 sxyz_tmp(3), uvsloc_tmp(2)
c
c  get bounding box
c
      bsize = 0.0d0
      nds = 12
      call get_bsize(nds, npts, srcvals, ndtarg, ntarg, targs, bsize)
c
c  intialize target info
c

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        flags(i) = -1
        sxyz(1:3,i) = 0.0d0
        ipatch_id_targ(i) = -1
        uvsloc(1:2,i) = 0
        dists(i) = 10*bsize
      enddo
C$OMP END PARALLEL DO

      allocate(cms(3,npatches), rads(npatches))
      
      call get_centroid_rads(npatches, norders, ixyzs, iptype, npts,
     1  srccoefs, cms, rads)
      do i=1,npatches
        rads(i) = rads(i)*rfac
      enddo
      
      nnz = 0
      call findnearmem(cms, npatches, rads, ndtarg, targs, ntarg, nnz)

      allocate(row_ptr(ntarg+1), col_ind(nnz))
      call findnear(cms, npatches, rads, ndtarg, targs, ntarg, row_ptr, 
     1  col_ind)
      
      
      
c
c  find the closest points
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, ipatch, rmin, ipts_init)
C$OMP$PRIVATE(l, rr, ntarg_tmp, flags_tmp, dist_tmp, sxyz_tmp)
C$OMP$PRIVATE(uvsloc_tmp)
      do i=1,ntarg
        dists(i) = 10*bsize 
        do j=row_ptr(i),row_ptr(i+1)-1
           ipatch = col_ind(j)
           rmin = 10*bsize 
           ipts_init = 1
           do l=ixyzs(ipatch),ixyzs(ipatch+1)-1
             rr = (targs(1,i) - srcvals(1,l))**2 + 
     1            (targs(2,i) - srcvals(2,l))**2 + 
     2            (targs(3,i) - srcvals(3,l))**2
             if(rr.le.rmin) then
               rmin = rr
               ipts_init = l
             endif
           enddo
           ntarg_tmp = 1
           flags_tmp = 1
           dist_tmp = 10*bsize
           call findnear_surface_point(ntarg_tmp, ndtarg, targs(1,i),
     1       npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals,
     2       ipatch_id, uvs_src, ipts_init, maxiter, tol, sxyz_tmp,
     3       uvsloc_tmp, dist_tmp, flags_tmp)
           if (dist_tmp.le.dists(i)) then
             dists(i) = dist_tmp
             sxyz(1:3,i) = sxyz_tmp(1:3)
             uvsloc(1:2,i) = uvsloc_tmp(1:2)
             flags(i) = flags_tmp
             ipatch_id_targ(i) = ipatch
           endif
        enddo
      enddo
C$OMP END PARALLEL DO      

      return
      end
    

      subroutine findnear_surface_point(ntargs, ndtarg, targs,
     1     npatches, norders, ixyzs, iptype, npts, srccoefs, 
     2     srcvals, ipatch_id, uvs_src, ipts_init, maxiter, tol, sxyz, 
     3     uvsloc, dists, flags)
c
c  Given a collection of targets, a surface, a patch corresponding
c  to each target, and an intial
c  guess for a closest discretization point on that patch, 
c  this subroutine identifies the closest (u,v) coordinate on that 
c  patch using Newton's method for minimizing d(x, \Gamma_{j})^2
c  where \Gamma_{j} is the patch of the initial guess 
c  for target x
c
c  Input arguments:
c    - ntargs: integer *8
c        number of targets
c    - ndtarg: integer *8
c        leading dimension of target array, should at least be 3.
c    - targs: real *8(ndtarg,ntargs)
c        The first three components must be the
c        xyz coordinates of the targets
c    - npatches: integer *8
c        number of patches
c    - norders: integer *8(npatches)
c        order of discretization on each patch
c    - ixyzs: integer *8(npatches+1)
c        ixyzs(i) denotes the starting location in srccoefs,
c        and srcvals array where information for patch i begins
c    - iptype: integer *8(npatches)
c        type of patch
c       * iptype = 1, triangular patch with RV nodes
c       * iptype = 11, quad patch with GL nodes
c       * iptype = 12, quad patch with Cheb nodes
c    - npts: integer *8
c        total number of discretization points on the boundary
c    - srccoefs: double precision (9,npts)
c        basis expansion coefficients of x, $\partial_{u} x$,
c        and $\partial_{v} x$.
c    - srcvals: double precision (12,npts)
c        x, $\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
c        discretization nodes
c    - ipatch_id: integer *8(npts)
c        patch id of the discretization points
c    - uvs_src: real *8(2,npts)
c        local uv coordinates of the discretization points
c    - ipts_init: integer *8(ntargs)
c        index of initial guess for point on surface. The distance
c        will be minimized on the patch that iptlocs lies on
c    - maxiter: integer *8
c        maximum number of iterations to run newton for
c    - tol: real *8
c        tolerance to be used for the Newton algorithm
c 
c  Output arguments:
c    - sxyz: real *8 (3,ntargs)
c        xyz coordinates of the closest point on the patch 
c        corresponding to the target-patch input
c    - uvsloc: real *8(2,ntargs)
c        uv coordinates of the closest point on the patch
c    - dists: real *8(ntargs)
c        distance of the closest point to the target
c    - flags: integer *8(ntargs)
c        flag corresponding to succesful execution of newton
c        flag(i) = 0 => Newton for target i converged successfully
c        otherwise Newton failed to achieve given tolerance for 
c        target i
c-------------------------------------

      
      implicit none
c     List of calling arguments
      integer *8, intent(in) :: ntargs, ndtarg
      real *8, intent(in) :: targs(ndtarg,ntargs)
      
      integer *8, intent(in) :: npatches, norders(npatches)
      integer *8, intent(in) :: ixyzs(npatches+1), iptype(npatches)
      integer *8, intent(in) :: npts
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      integer *8, intent(in) :: ipatch_id(npts)
      real *8, intent(in) :: uvs_src(2,npts)

      integer *8, intent(in) :: ipts_init(ntargs)

      
      integer *8, intent(in) :: maxiter
      real *8, intent(in) :: tol

      real *8, intent(out) :: sxyz(3,ntargs), uvsloc(2,ntargs)
      real *8, intent(out) :: dists(ntargs)
      integer *8, intent(out) :: flags(ntargs)
      
c     List of local variables
      integer *8 nomax, npmax
      real *8, allocatable :: rsc1(:,:), rat1(:,:), rat2(:,:)
      integer *8, allocatable :: inorderuse(:),iptypeuse(:)
      integer *8, allocatable :: iuni(:,:), iuniind(:)

      
      integer *8 snpols,snpols2,nn,it,lpatchidx
      integer *8 ipnt, idx, k, i, iffast, iind, ipatch
      integer *8 ipt, ipuse, istart, itarg, nouse, npols
      integer *8 nuni
      real *8 sxyz0(3),bs,uv0(2)
      real *8, allocatable :: pols(:)

      nomax = maxval(norders)
      npmax = (nomax+1)*(nomax+1)
      allocate(pols(npmax))

c
c  identify unique pairs of norder, iptype combos
c

      allocate(inorderuse(ntargs),iptypeuse(ntargs))
      allocate(iuni(2,ntargs),iuniind(ntargs))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ipatch)
      do i=1,ntargs
        ipatch = ipatch_id(ipts_init(i))

        inorderuse(i) = norders(ipatch)
        iptypeuse(i) = iptype(ipatch) 
        iuni(1,i) = 0
        iuni(2,i) = 0
        iuniind(i) = 0
      enddo
C$OMP END PARALLEL DO
      nuni = 0
      call get_iuni2(ntargs, inorderuse, iptypeuse, nuni, iuni, 
     1  iuniind)
      allocate(rsc1(npmax,nuni), rat1(2*(nomax+1),nuni))
      allocate(rat2(3*npmax,nuni))


      do i=1,nuni
        nouse = iuni(1,i)
        ipuse = iuni(2,i)

        if(ipuse.eq.1) then
          call koornf_init(nouse, rat1(1,i), rat2(1,i), rsc1(1,i))
        endif
      enddo

      iffast = 1
      
c     loop over target points
      do itarg = 1,ntargs
         ipt = ipts_init(itarg)
c
c  get initial uv coordinates and patch id 
c
         uv0(1) = uvs_src(1,ipt)
         uv0(2) = uvs_src(2,ipt)
         ipatch = ipatch_id(ipt)
         istart = ixyzs(ipatch)
         npols = ixyzs(ipatch+1) - ixyzs(ipatch)
         iind = iuniind(itarg)
         
c     Newton to find closest point, returns point (u,v) in parameter space
         call newton_nearp(norders(ipatch), npols, iptype(ipatch),
     1        srccoefs(1,istart), srcvals(1,istart), targs(1,itarg),
     2        uv0, tol, maxiter, uvsloc(1,itarg), flags(itarg), iffast,
     3        rat1(1,iind), rat2(1,iind), rsc1(1,iind))
         
c     Convert (u,v) in parameter space to a point sxyz in R^3
         call get_basis_pols(uvsloc(1,itarg), norders(ipatch), npols,
     1     iptype(ipatch), pols)
         sxyz(1,itarg) = 0.0d0
         sxyz(2,itarg) = 0.0d0
         sxyz(3,itarg) = 0.0d0        
         do k=1,npols
            sxyz(1,itarg) = sxyz(1,itarg) + 
     1          srccoefs(1,istart+k-1)*pols(k)
            sxyz(2,itarg) = sxyz(2,itarg) + 
     1          srccoefs(2,istart+k-1)*pols(k)
            sxyz(3,itarg) = sxyz(3,itarg) + 
     1          srccoefs(3,istart+k-1)*pols(k)
         end do
         
c     Compute distance between target point its closest point sxyz
         dists(itarg)=(sxyz(1,itarg)-targs(1,itarg))**2
     1        + (sxyz(2,itarg)-targs(2,itarg))**2
     2        + (sxyz(3,itarg)-targs(3,itarg))**2
         dists(itarg) = dsqrt(dists(itarg))
      enddo

      return    
      end 
c     
c     
c      
c
      subroutine newton_nearp(norder, npols, iptype, srccoefs, srcvals, 
     1     pt, uv0, tol, maxiter, uvs, flag, iffast, rat1, rat2, rsc1)
c
c  This subroutine runs newton to find the closest point on a patch.
c  
c  Input arguments:
c    - norder: integer *8
c        order of patch
c    - npols: integer *8
c        number of discretization nodes on the patch
c    - iptype: integer *8(npatches)
c        type of patch
c       * iptype = 1, triangular patch with RV nodes
c       * iptype = 11, quad patch with GL nodes
c       * iptype = 12, quad patch with Cheb nodes
c    - srccoefs: double precision (9,npols)
c        basis expansion coefficients of x, $\partial_{u} x$,
c        and $\partial_{v} x$.
c    - srcvals: double precision (12,npols)
c        x, $\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
c        discretization nodes
c    - pt: real *8(3)
c        xyz coordinates of the target
c    - uv0: real *8(2)
c        uv coordinates of the initial guess for newton (Note this
c        variable is overwritten in this routine and is equal 
c        to the output variable uvs at the end of it)
c    - tol: real *8
c        precision to run newton to
c    - maxiter: integer *8
c        maximum number of newton iterations
c    - iffast: integer *8
c        flag for using precomputed coefficients for basis function
c        evaluation
c    - rat1: real *8 (*)
c        recurrence coefficients 1
c    - rat2: real *8 (*)
c        recurrence coefficients 2
c    - rsc1: real *8 (*)
c        scaling factors
c 
c  Output arguments:
c    - uvs: real *8(2)
c        uv coordinates of the point on the patch which minimizes the
c        distance
c    - flags: integer *8(ntargs)
c        flag corresponding to succesful execution of newton
c        flag = 0 => Newton for target converged successfully
c        otherwise Newton failed to achieve given tolerance 


      implicit none

!List of calling arguments
      integer *8, intent(in) :: norder
      integer *8, intent(in) :: npols
      integer *8, intent(in) :: iptype
      integer *8, intent(in) :: maxiter
      real ( kind = 8 ), intent(in) :: srcvals(12,npols)
      real ( kind = 8 ), intent(in) :: srccoefs(9,npols)
      real ( kind = 8 ), intent(in) :: pt(3)
      real ( kind = 8 ), intent(in) :: tol
      real ( kind = 8 ), intent(inout) :: uv0(2)
      integer *8, intent(in) :: iffast
      real ( kind = 8), intent(in) :: rat1(*)
      real ( kind = 8), intent(in) :: rat2(*)
      real ( kind = 8), intent(in) :: rsc1(*)
      real ( kind = 8 ), intent(out) :: uvs(2)
      integer *8, intent(out) :: flag
      

!List of local variables
      real ( kind = 8 ) err, fval, fval_new, alpha, descent
      integer *8 iter, ibt, islarge
      real ( kind = 8 ) gradf(2), hesf(2,2), delta(2), dethesf
      real ( kind = 8 ) uvtmp(2)
      real *8, allocatable :: ptmp(:), dtmp(:,:)

      allocate(ptmp(npols), dtmp(2,npols))

      err = (tol+1.0d0)**2
      iter = 0
      uvs(1) = uv0(1)
      uvs(2) = uv0(2)
      flag = 1

      do while ((iter < maxiter) .and. flag.eq.1)

         call eval_fgradhess(norder, npols, iptype, srccoefs, srcvals,
     1        uv0, pt, fval, gradf, hesf, ptmp, dtmp, iffast, rat1,
     2        rat2, rsc1)

         dethesf = hesf(1,1)*hesf(2,2) - hesf(1,2)*hesf(2,1)
         delta(1) = (-hesf(1,2)*gradf(2) + hesf(2,2)*gradf(1))/dethesf
         delta(2) =  (hesf(1,1)*gradf(2) - hesf(2,1)*gradf(1))/dethesf

         descent = delta(1)*gradf(1) + delta(2)*gradf(2)
         if (descent.lt.0) then
           delta(1) = gradf(1)
           delta(2) = gradf(2)
           descent = gradf(1)**2 + gradf(2)**2
         endif

         uvtmp(1:2) = uv0(1:2) - delta(1:2)

c        Backtracking line search (Armijo condition)
         alpha = 1.0d0
         do ibt = 1, 30
            uvtmp(1) = uv0(1) - alpha*delta(1)
            uvtmp(2) = uv0(2) - alpha*delta(2)
            call eval_fgradhess(norder, npols, iptype, srccoefs, 
     1         srcvals, uvtmp, pt, fval_new, gradf, hesf, ptmp, dtmp,
     2         iffast, rat1, rat2, rsc1)

            if (fval_new <= fval - 0.3d0 * alpha * descent) exit
            alpha = 0.2d0 * alpha
         enddo

c
c  project back
c
         if (iptype.eq.1) then
           if (uvtmp(1)+uvtmp(2).gt.1.0d0) 
     1         uvtmp(1:2) = uvtmp(1:2)/(uvtmp(1) + uvtmp(2))
           if (uvtmp(1).lt.0) uvtmp(1) = 0
           if (uvtmp(2).lt.0) uvtmp(2) = 0
         else
           if (uvtmp(1).lt.-1) uvtmp(1) = -1.0d0
           if (uvtmp(1).gt.1) uvtmp(1) = 1.0d0
           if (uvtmp(2).lt.-1) uvtmp(2) = -1.0d0
           if (uvtmp(2).gt.1) uvtmp(2) = 1.0d0
         endif


         uvs(1) = uvtmp(1)
         uvs(2) = uvtmp(2)


         err = (uvs(1)-uv0(1))**2 + (uvs(2)-uv0(2))**2

         uv0(1) = uvs(1)
         uv0(2) = uvs(2)
         if (err.le.tol**2. or. descent.le.tol**2) then
           flag = 0
         else
           flag = 1
         endif
         iter = iter + 1
      enddo

      return
      end subroutine newton_nearp
!     
!     
!     
!     
!     
      subroutine eval_fgradhess(norder, npols, iptype, srccoefs,
     1  srcvals, uvs, pt, fval, gradf, hesf, pols, ders, iffast, rat1,
     2  rat2, rsc1)
      implicit none

!List of calling arguments
      integer *8, intent(in) :: norder
      integer *8, intent(in) :: npols
      integer *8, intent(in) :: iptype
      real ( kind = 8 ), intent(in) :: srcvals(12,npols)
      real ( kind = 8 ), intent(in) :: srccoefs(9,npols)
      real ( kind = 8 ), intent(in) :: pt(3)
      real ( kind = 8 ), intent(in) :: uvs(2)
      integer *8, intent(in) :: iffast
      real ( kind = 8 ), intent(in) :: rat1(*)
      real ( kind = 8 ), intent(in) :: rat2(*)
      real ( kind = 8 ), intent(in) :: rsc1(*)
      real ( kind = 8 ), intent(out) :: fval
      real ( kind = 8 ), intent(out) :: gradf(2)
      real ( kind = 8 ), intent(out) :: hesf(2,2)
      real ( kind = 8 ) :: pols(npols), ders(2,npols)

!List of local variables
      integer *8 i
      real ( kind = 8 ) x,y,z,xu,yu,zu,xv,yv,zv
      real ( kind = 8 ) xuu,xuv,xvv,yuu,yuv,yvv,zuu,zuv,zvv,aux1(2)


      if (iffast.ne.1) then
        call get_basis_ders(uvs, norder, npols, iptype, pols, ders)
      else
        if (iptype.eq.1) then
          call koornf_ders(uvs, norder, npols, pols, ders, rat1, rat2, 
     1      rsc1)
        else
          call get_basis_ders(uvs, norder, npols, iptype, pols, ders)
        endif
      endif
      x = 0
      y = 0
      z = 0
      xu = 0
      xv = 0
      xuu = 0
      xuv = 0
      xvv = 0
      yu = 0
      yv = 0
      yuu = 0
      yuv = 0
      yvv = 0
      zu = 0
      zv = 0
      zuu = 0
      zuv = 0
      zvv = 0

      do i=1,npols
         x = x + srccoefs(1,i)*pols(i)
         y = y + srccoefs(2,i)*pols(i)
         z = z + srccoefs(3,i)*pols(i)

         xu = xu + srccoefs(4,i)*pols(i)
         yu = yu + srccoefs(5,i)*pols(i)
         zu = zu + srccoefs(6,i)*pols(i)

         xv = xv + srccoefs(7,i)*pols(i)
         yv = yv + srccoefs(8,i)*pols(i)
         zv = zv + srccoefs(9,i)*pols(i)

         xuu = xuu + srccoefs(4,i)*ders(1,i)
         yuu = yuu + srccoefs(5,i)*ders(1,i)
         zuu = zuu + srccoefs(6,i)*ders(1,i)

         xuv = xuv + srccoefs(4,i)*ders(2,i)
         yuv = yuv + srccoefs(5,i)*ders(2,i)
         zuv = zuv + srccoefs(6,i)*ders(2,i)

         xvv = xvv + srccoefs(7,i)*ders(2,i)
         yvv = yvv + srccoefs(8,i)*ders(2,i)
         zvv = zvv + srccoefs(9,i)*ders(2,i)
      enddo


      fval = (x-pt(1))**2 + (y-pt(2))**2 + (z-pt(3))**2
      gradf(1)=2*((x-pt(1))*xu+(y-pt(2))*yu+(z-pt(3))*zu)
      gradf(2)=2*((x-pt(1))*xv+(y-pt(2))*yv+(z-pt(3))*zv)
      hesf(1,1)=2*(xu**2+yu**2+zu**2)+2*((x-pt(1))*xuu+(y-pt(2))*yuu
     1     +(z-pt(3))*zuu)
      hesf(1,2)=2*(xu*xv+yu*yv+zu*zv)+2*((x-pt(1))*xuv+(y-pt(2))*yuv
     1     +(z-pt(3))*zuv)
      hesf(2,1)=hesf(1,2)
      hesf(2,2)=2*(xv**2+yv**2+zv**2)+2*((x-pt(1))*xvv+(y-pt(2))*yvv
     1     +(z-pt(3))*zvv)


      return
      end subroutine eval_fgradhess
!     
!     
!     
!     
!     
!     
      
