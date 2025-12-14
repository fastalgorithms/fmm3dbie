!
!  This file contains codes for extracting precomputed special-
!  purpose quadratures for evaluating integrals of the form
!       \int_{S} K(x,y) \sigma(y) dS_{y}           --- (1)
!
!  where S is a convex polygon
!
!  The code uses quadratures generated in
!  https://github.com/zgimbutas/selfquad3d.git available 
!  under the public domain license
!
!  The code is an adaptation of the routines in selfquad3d.f
!  in the above repository
!
! 
      subroutine self_quadrature(norder, ipv, verts, nv, x0, y0, &
         dr, nquad, xs, ys, whts, ier)
!
!  This routine extracts pregenerated generalized quadrature rules
!  for integrating functions of the form
!
!     \int_{S} K(x,y) \sigma(y) dS_{y}
!
!  If ipv = 0, then 
!       K(x,y) = log(|x-y|) smooth + 1/|x-y| smooth + smooth
!  If ipv = 1
!       K(x,y) = 1/|x-y| smooth + p.v. 1/|x-y|^2*smooth
!  If ipv = 2
!       K(x,y) = 1/|x-y| smooth + p.v. 1/|x-y|^2*smooth + 
!         f.p. 1/|x-y|^3*smooth (under construction)
!
!  Here S is a convex polygon in R^2 whose boundary is defined
!  by an ordered set of vertices
!
!  Input arguments:
!    - norder: integer
!        order of polynomials for representing sigma on S
!    - ipv: integer
!        nature of singularity of the kernel K (see above)
!        * ipv = 0, single layer like + log singular kernels
!        * ipv = 1, p.v. integrals
!        * ipv = 2, hypersingular integrals
!    - verts: real *8(2,nv)
!        coordinates defining the convex polygon
!    - nv: integer
!        number of vertices
!    - x0: real *8
!        x coordinate of target point on the polygon S
!    - y0: real *8
!        y coordiante of target point on the polygon S
!    - dr: real *8(3,2)
!        u and v derivatives of the surface parametrization
!        at \rho(x0, y0) where \rho: S\to \Gamma
!        is a parametrization of the surface
!
!  Output arguments:
!    - nquad: integer
!        number of nodes in the final quadrature rule
!    - xs: real *8(nquad)
!        x coordinates of quadrature nodes on S
!    - ys: real *8(nquad)
!        y coordinates of quadrature nodes on S
!    - whts: real *8(nquad)
!        quadrature weights on S
!    - ier: integer
!        error code
!        ier = 0 => successful execution
!
      implicit none
      integer, intent(in) :: norder
      integer, intent(in) :: ipv
      integer, intent(in) :: nv
      real *8, intent(in) :: verts(2,nv), x0, y0
      real *8, intent(in) :: dr(3,2)
      integer, intent(out) :: nquad, ier
      real *8, intent(out) :: xs(*), ys(*), whts(*)

!
!  temporary variables
!
      real *8 qdr(3,2), rdr(2,2), drdetinv
      real *8 vtransf(2,nv+1), projs(2,nv), thetas(nv)
      real *8 hdis(nv)
      real *8 d, rtmp1, rtmp2, vtmp(2)
      real *8 r, w, rcut, rcmin, t, tw, dd 
      integer i, j, k, l, nq1, nlege, ifwhts, nthet, iquad 
      real *8 xr(150), wr(150), pi, xtmp, ytmp
      integer nr, nr0
      integer ier0
      data pi/3.14159265358979323846264338327950288d0/

!
!  compute qr decomposition of jacobian vector.
!  The jacobian tends to be reasonably conditioned, i.e. 
!  better than 1/sqrt(machine precision), so classical
!  Gram-Schmidt suffices
!
      qdr(1:3,1:2) = dr(1:3,1:2)
      rdr(1,1) = sqrt(qdr(1,1)**2 + qdr(2,1)**2 + qdr(3,1)**2)
      qdr(1:3,1) = qdr(1:3,1)/rdr(1,1)

      rdr(1,2) = qdr(1,1)*qdr(1,2) + qdr(2,1)*qdr(2,2) + &
        qdr(3,1)*qdr(3,2)
      qdr(1:3,2) = qdr(1:3,2) - rdr(1,2)*qdr(1:3,1)

      rdr(2,2) = sqrt(qdr(1,2)**2 + qdr(2,2)**2 + qdr(3,2)**2)
      qdr(1:3,2) = qdr(1:3,2)/rdr(2,2)

      drdetinv = 1.0d0/abs(rdr(1,1)*rdr(2,2))

!
!  compute vertices in conformal coordinates, i.e.
!  in basis vectors given by q
!

      do i=1,nv
        rtmp1 = verts(1,i) - x0
        rtmp2 = verts(2,i) - y0
        
        vtransf(1,i) = rdr(1,1)*rtmp1 + rdr(1,2)*rtmp2
        vtransf(2,i) = rdr(2,2)*rtmp2
      enddo
      vtransf(1:2,nv+1) = vtransf(1:2,1)

!
!  find projections to the sides and the cut-out radius
!  for principal value and hypersingular quadratures
!
      vtmp(1) = 0
      vtmp(2) = 0
      do i=1,nv
        call comp_projs_selfquad(vtmp, vtransf(1,i), vtransf(1,i+1), &
          hdis(i), projs(1,i), thetas(i))
      enddo
      rcmin = minval(hdis)
      rcut = rcmin/2

!
!  Now start assembling on the triangles. Ignore almost collinear triangles,
!  for compact quadratures, this can then inturn be used to evaluate
!  layer potentials on the boundary of the triangle
!  
!
      nquad = 0
      ier = 0
      nq1 = 0
      do i = 1,nv
        if(abs(thetas(i)).gt.3.11d0) then
          call get_ggq_triangle_origin_quadrature(norder, ipv, rcut, &
            vtransf(1,i), projs(1,i), xs(1+nquad), ys(1+nquad), &
            whts(1+nquad), nq1, ier0)
          nquad = nquad + nq1
          ier = max(ier, ier0)

          call get_ggq_triangle_origin_quadrature(norder, ipv, rcut, &
            projs(1,i), vtransf(1,i+1), xs(1+nquad), ys(1+nquad), &
            whts(1+nquad), nq1, ier0)
          nquad = nquad + nq1
          ier = max(ier, ier0)
        else
          call get_ggq_triangle_origin_quadrature(norder, ipv, rcut, &
            vtransf(1,i), vtransf(1,i+1), xs(1+nquad), ys(1+nquad), &
            whts(1+nquad), nq1, ier0)
          nquad = nquad + nq1
          ier = max(ier, ier0)
        endif
      enddo
!
!  Add in tensor product quadrature on the cut out disk
!
!
      nthet = 3*(norder+1) + 2
      nlege = (norder/2) + 1
      if(ipv.eq.1) then
        ifwhts = 1
        call legewhts(nlege, xr, wr, ifwhts)
        nr = nlege
      elseif (ipv.eq.2) then
        call hs_disk_quad(nlege, xr, wr, nr)
      endif
      
      if (ipv.eq.1.or.ipv.eq.2) then
        dd = (2*pi)/nthet
        do i = 1,nr
          r = rcut/2*xr(i) + rcut/2
          w = rcut/2*wr(i)
          do j = 1,nthet
            iquad = nquad + (i-1)*nthet + j
            t = dd*(j-1)
            tw = dd
            xs(iquad) = r*cos(t)
            ys(iquad) = r*sin(t)
            whts(iquad) = tw*r*w
          enddo
        enddo
        nquad = nquad + nthet*nr
      endif
!
!  Inverse qr mapping to all the quadrature nodes
! 
!
      do i = 1,nquad
        ytmp = ys(i)/rdr(2,2)
        xtmp = xs(i) - rdr(1,2)*ytmp
        xtmp = xtmp/rdr(1,1)
        
        xs(i) = xtmp + x0
        ys(i) = ytmp + y0
        whts(i) = whts(i)*drdetinv
      enddo


      ier = 0
      return
      end subroutine

!
!
!
!
!
      subroutine comp_projs_selfquad(v0, v1, v2, h, proj, thet)
!
!  Find the signed length of height and projections of v0
!  to the line segment defined by v1 and v2
!
!  Input arguments:
!    - v0: real *8(2)
!        vertex whose projection is to be computed
!    - v1: real *8(2)
!        start vertex of line segment
!    - v2: real *8(2)
!        end vertex of line segment
!
!  Output arguments:
!    - h: real *8
!        signed distance of v0, to the line segment formed
!        by v1 -- v2
!    - proj: real *8(2)
!        location on line segment formed by v1 -- v2, such
!        that v0 - proj \perp v1 - v2
!    - thet: real *8
!        angle formed at v0, by the triangle v1 -- v0 -- v2
!
      implicit none
      real *8, intent(in) :: v0(2), v1(2), v2(2)
      real *8, intent(out) :: h, proj(2), thet

!
!  Temporary variables
!
      real *8 v01(2), dir(2), coeff, rn(2), dvec(2)
      real *8 v20(2), dx, dy, dsinv
      complex *16 z10bar, z20, ztmp, ima
      data ima/(0.0d0,1.0d0)/
      
      
      v01(1) = v0(1) - v1(1)
      v01(2) = v0(2) - v1(2)

      z10bar = v1(1) - v0(1) - ima*(v1(2) - v0(2))
      

      v20(1) = v2(1) - v0(1)
      v20(2) = v2(2) - v0(2)

      z20 = v2(1) - v0(1) + ima*(v2(2) - v0(2))

      ztmp = z10bar*z20
      thet = atan2(imag(ztmp), real(ztmp))

      dir(1) = v2(1) - v1(1)
      dir(2) = v2(2) - v1(2)

      rn(1) = -dir(2)
      rn(2) = dir(1)

      dsinv = 1.0d0/sqrt(dir(1)**2 + dir(2)**2)

      coeff = (v01(1)*dir(1) + v01(2)*dir(2))*dsinv**2

      proj(1:2) = v1(1:2) + coeff*dir(1:2)
      dvec(1:2) = v0(1:2) - proj(1:2)
      h = (dvec(1)*rn(1) + dvec(2)*rn(2))*dsinv

      return
      end subroutine
!
!
!
!
!
      subroutine get_ggq_triangle_origin_quadrature(norder, ipv, rcut, &
        v0, v1, xs, ys, ws, nquad, ier)
!
!  Return a quadrature for evaluating radially singular functions of
!  the form
!     \int_{T} K(x,y) \sigma(y) dT
!
!  If ipv = 0, then 
!       K(x,y) = log(|x-y|) smooth + 1/|x-y| smooth + smooth
!  If ipv = 1
!       K(x,y) = 1/|x-y| smooth + p.v. 1/|x-y|^2*smooth
!  If ipv = 2
!       K(x,y) = 1/|x-y| smooth + p.v. 1/|x-y|^2*smooth + 
!         f.p. 1/|x-y|^3*smooth (under construction)
!
!  Here T is a triangle with vertices (0,0), v0, v1.
!
!  Input arguments:
!    - norder: integer
!        order of polynomials for representing sigma on S
!    - ipv: integer
!        nature of singularity of the kernel K (see above)
!        * ipv = 0, single layer like + log singular kernels
!        * ipv = 1, p.v. integrals
!    - rcut: real *8
!        radius of disk around the origin for which the quadrature
!        is constructed separately (only applicable for ipv>=1)
!    - v0: real *8(2)
!        vertex of the triangle
!    - v1: real *8(2)
!        the other vertex of the triangle
!
!  Output arguments:
!    - xs: real *8(nquad)
!        x coordinates of quadrature nodes on the triangle
!    - ys: real *8(nquad)
!        y coordinates of quadrature nodes on the triangle
!    - ws: real *8(nquad)
!        quadrature weights
!    - ier: integer
!        error code
!        * ier = 0 => successful execution
!        * ier = 4 => couldn't fetch theta quadratures
!
!
      implicit none
      integer, intent(in) :: norder, ipv
      real *8, intent(in) :: rcut, v0(2), v1(2)
      integer, intent(out) :: nquad, ier
      real *8, intent(out) :: xs(*), ys(*), ws(*)
!
! Temporary variables
!
      real *8 r0, r1, xmeas, rr, sina, cosa, xx, yy
      real *8 rcut_scaled, r, theta
      integer iquad, ier0
      integer ifreflect, irad, nr, nr0, nt, i, j, it, istart, ir
      real *8 xr(200), wr(200), xt(200), wt(200)
      real *8 xtmp, ytmp
      real *8 ct, st, t, twht, rwht, ruse, rsc, rr2, rl
      
      nquad = 0
      ier = 0

!
! Test for collinearity
!
      xmeas = 1.0d0
      if ((abs(v0(1)).gt.abs(v0(2))).and.(abs(v1(1)).gt.abs(v1(2)))) then
        xmeas = abs(v0(2)/v0(1) - v1(2)/v1(1))
      else
        xmeas = abs(v0(1)/v0(2) - v1(1)/v1(2))
      endif

      if (xmeas.le.1.0d-8) return

!
!  Points are no longer collinear
!
      r0 = sqrt(v0(1)**2 + v0(2)**2)
      r1 = sqrt(v1(1)**2 + v1(2)**2)
!
!  Place the longer edge on the y=0 axis, rotate and
!  reflect about y=0 axis if needed
!
      if (r0.gt.r1) then
        rr = r0
        sina = v0(2)/r0
        cosa = v0(1)/r0
        xx = cosa*v1(1) + sina*v1(2)
        yy = -sina*v1(1) + cosa*v1(2)
      else
        rr = r1
        sina = v1(2)/r1
        cosa = v1(1)/r1
        xx = cosa*v0(1) + sina*v0(2)
        yy = -sina*v0(1) + cosa*v0(2)
      endif

      ifreflect = 0
      if (yy.lt.0) then
        yy = -yy
        ifreflect = 1
      endif
!
!  normalize
!
      xx = xx/rr
      yy = yy/rr
      rcut_scaled = rcut/rr
      if (ipv.eq.0) rcut_scaled = 0

      r = sqrt(xx**2 + yy**2)
      theta = atan2(yy, xx)

!
!  get which theta quadratures to fetch
!  determine irad, based on norder
!
      if (norder.lt.2) then
        irad = 1
      elseif (norder.lt.4.and.norder.ge.2) then
        irad = 2
      elseif (norder.lt.6.and.norder.ge.4) then
        irad = 3
      elseif (norder.lt.8.and.norder.ge.6) then
        irad = 4
      elseif (norder.lt.10.and.norder.ge.8) then
        irad = 5
      elseif (norder.lt.12.and.norder.ge.10) then
        irad = 6
      elseif (norder.lt.16.and.norder.ge.12) then
        irad = 7
      elseif (norder.lt.20.and.norder.ge.16) then
        irad = 8
      else
        ier = 4
        print *, "norder too high for theta quadratures"
        print *, "Returning without returning self quadrature"
        return
      endif


!
!
      ir = 0
      it = 0
      call get_radfetch_param(ipv, irad, r, theta, nt, ir, it, ier)
      if (ier.gt.0) then
        print *, "Couldn't find theta quadratures"
        print *, "Returning without returning self quadrature"
        return
      endif
      call radfetch(ipv, irad, ir, it, nt, xt, wt)
!
!  If ipv = 0, then the radial quadratures are indepdent of
!  the transverse quadratures. Pregenerate them here
!
      if(ipv.eq.0) then
        nr0 = norder + 4
        call rquads_rlogr_rinv(nr0, xr, wr, nr)
      endif

      istart = 0
      do j = 1,nt
        t = xt(j)*theta
        twht = wt(j)*theta
        ct = cos(t)
        st = sin(t)
          
        ruse = r*sin(theta)/(r*sin(theta-t) + st) 
        rl = ruse - rcut_scaled

        if (ipv.eq.1.or.ipv.eq.2) then
          call rquads_pv_hs(norder, ipv, rcut_scaled, rl, xr, wr, nr, &
              ier)
          if (ier.gt.0) then
            print *, "Failed to extract radial quadratures"
            print *, "Returning without returning self quadrature"
            return
          endif
        endif
!
!  Now construct the tensor product quadrature
! 
        do i = 1,nr
          iquad = istart + i
          rsc = rcut_scaled + xr(i)*rl
          rwht = wr(i)*rl
          xs(iquad) = rsc*ct
          ys(iquad) = rsc*st
          ws(iquad) = rsc*rwht*twht
        enddo
        istart = istart + nr
      enddo

      nquad = istart 
!
!  Now undo the scaling and the transformations at the top
!
      rr2 = rr*rr
      do i = 1,nquad
        xs(i) = xs(i)*rr
        ys(i) = ys(i)*rr
        ws(i) = ws(i)*rr2

        if (ifreflect .eq. 1) ys(i) = -ys(i)
        xtmp = cosa*xs(i) - sina*ys(i)
        ytmp = sina*xs(i) + cosa*ys(i)
        xs(i) = xtmp
        ys(i) = ytmp
      enddo

      return
      end subroutine



      subroutine rquads_rlogr_rinv(norder, xs, ws, nr)
      implicit none
      integer norder, nr
      real *8 xs(*), ws(*)

      if (norder.gt.24) then
        print *, "no quadrautres found in rquads_rlogr_rinv"
        print *, "returning without any quadratures"
        nr = 0
        return
      endif
      
      INCLUDE 'ggq-self-quads/rlogr_rinv_quads.txt'
      nr = norder

      return
      end


      subroutine rquads_pv_hs(norder, ipv, a, b, xr, wr, nr, ier)
!
!  This subroutine returns special-purpose quadratures for the integrals
!  \int_{0}^{b} f(x) dx where   
!  
!  f(x) = phi(x) + 1/(x+a)    when ipv = 1, and
!  f(x) = phi(x) + 1/(x+a) + 1/(x+a)^2 when ipv = 2 
!
!  Input arguments:
!    - norder: integer
!        order of polynomials for representing sigma on S
!    - ipv: integer
!        nature of singularity of the kernel K (see above)
!        * ipv = 1, p.v. integrals
!        * ipv = 2, hypersingular integrals
! 
!  Output arguments:
!    - xr: real *8(nr)
!        quadrature nodes
!    - wr: real *8(nr)
!        quadarture weights
!    - nr: integer
!        number of quadrature nodes/weights
!    - ier: integer
!        error code
!        ier = 0 => successful execution
!----------------------------------
      implicit none
      integer, intent(in) :: ipv, norder
      real *8, intent(in) :: a, b
      real *8, intent(out) :: xr(*), wr(*)
      integer, intent(out) :: nr, ier

!
!  temporary variables
! 
      real *8 rs(29), delta
      integer i, n

      delta = a/b

      do i=1,20
        rs(i) = 10**(-20+i-1)
      enddo

      do i=21,29
        rs(i) = (i-19)*0.1d0
      enddo

      n = 0
      do i = 1,28
        if( delta .ge. rs(i) .and. delta .le. rs(i+1) ) n = i
      enddo

      if (n.eq.0) then
        nr = 0
        ier = 4
        return
      endif

      if (norder.gt.20) then
        ier = 4
        print *, "Order too high for self-quadrature"
        print *, "Returning without returning self-quadrature"
        return
      endif

      if (ipv.eq.1) then
        INCLUDE 'ggq-self-quads/rpv_quads.txt'
      endif

      if (ipv.eq.2) then
        INCLUDE 'ggq-self-quads/rhs_quads.txt'
      endif

      return
      end


      subroutine hs_disk_quad(n, xr, wr, nr)
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: n
      real *8, intent(out) :: xr(*), wr(*)
      integer, intent(out) :: nr
      
      INCLUDE 'ggq-self-quads/rhs_quads_disk.txt'

      return
      end
!
!
!
!
      subroutine get_radfetch_param(ipv, irad, r, t, n, ir, it, ier)
!---------------------
!
!  this subroutine determines the quadrature 
!  region number corresponding to parameters
!  r,t,and irad
!
!  Input arugments:
!
!    - ipv: integer
!        flag for extracting compact/pv quadratures
!        ipv = 0, compact quadratures
!        ipv = 1, pv quadratures
!    - irad: integer
!        type of radial quadrature to be used
!        for norder<=4, use irad=1,
!        for 4<norder<=8, use irad=2,
!        else use irad=3, there will be some loss in accuracy
!        if norder>12. Support for higher orders will be made
!        available with irad=4, and so on
!
!    - r: double precision
!         radial coordinate in triangle
!    - t: double precision
!         angle coordinate in triangle
!
!  Output arguments:
!
!    - n: integer
!        number of quadrature nodes in the region
!    - ir: integer
!        radial region number
!    - it: integer
!        theta region number
!
      implicit none
      integer, intent(in) :: irad,ipv,ier
      real *8, intent(in) :: r,t
      integer, intent(out) :: n,ir,it

      real *8 rs(2,25),ts(2,25),eps
      integer nn(25,25),nrs,nts
      integer i

      eps = 1.0d-14

      if(irad.lt.1.or.irad.gt.9) then
        call prinf('invalid parameter irad*',i,0)
      endif

      INCLUDE 'ggq-self-quads/radfetch-param-data.txt'

      do ir=1,nrs
        if (rs(1,ir) .le. r .AND. r .le. rs(2,ir)+eps) goto 1100
      enddo
      call prinf('r in invalid range in radfetch0*',i,0)
 1100 continue

      do it=1,nts
        if (ts(1,it) .le. t .AND. t .le. ts(2,it)+eps) goto 1200
      enddo
      call prinf('t in invalid range in radfetch0*',i,0)
 1200 continue

      n = nn(it,ir)

      return
      end
!
!
!
!
!
      
      subroutine radfetch(ipv, irad, ir, it, nt, ts, wts)
!---------------------
!
!  this subroutine extracts the quadrature 
!  corresponding to region number ipv,irad,ir,it
!
!  Input arugments:
!
!    - ipv: integer
!        flag for extracting compact/pv quadratures
!        ipv = 0, compact quadratures
!        ipv = 1, pv quadratures
!    - irad: integer
!        type of radial quadrature to be used
!        for norder<=4, use irad=1,
!        for 4<norder<=8, use irad=2,
!        else use irad=3, there will be some loss in accuracy
!        if norder>12. Support for higher orders will be made
!        available with irad=4, and so on
!    - ir: integer
!        radial region number
!    - it: integer
!        theta region number
!    - nt: integer
!        number of quadrature nodes in the region
!
!  Output arguments:
!    
!    - ts: double precision(nt)
!        quadrature nodes
!    - wts: double precision(nt)
!        quadrature weights
!
      implicit none
      integer, intent(in) :: irad,ir,it,nt,ipv
      real *8, intent(out) :: ts(nt),wts(nt)

      INCLUDE 'ggq-self-quads/radfetch-data.txt'

      return
      end
