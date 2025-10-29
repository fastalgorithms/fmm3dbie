      subroutine mesh_circle_pts_mem(rmid, npars, iort, norder, &
        iptype0, npatches, npts)
!
!  This subroutine returns the numbers of patches and points on
!  the mesh of a unit circle
!
!  There is middle square of length rmid, and four sectors going
!  from the vertices of the square to the boundary of the
!  circle in a radial fashion.
!
!  The middle square is discretized into (npars(1)*npars(1)) 
!  smaller squares, and each of the sectors is discretized 
!  into further (npars(2)*npars(3)) squares, with npars(2)
!  intervals in the radial direction and npars(3) 
!  intervals in the angular direction
!  
!
!  Input Arguments:
!    - rmid: real *8
!        radius of the inner square, should be between 0 and 1/sqrt(2)
!    - npars: integer(3)
!        * npars(1) = nmid, number of intervals in the middle square
!        * npars(2) = nr, number of intervals in the radial direction
!                     in each of the sectors
!        * npars(3) = nt, number of intervals in the angular direction
!                     in each of the sectors
!    - iort: integer
!        orientiation flag
!        if iort = 1, then mesh has normals pointed in +z direction
!        if iort = -1, then mesh has normals pointed in -z direction 
!    - norder: integer
!        order of discretization on each patch
!    - iptype0: integer
!        type of patch to be used in the discretization
!        * iptype = 1, triangular patch discretized using RV nodes
!        * iptype = 11, quadrangular patch discretized with GL nodes
!        * iptype = 12, quadrangular patch discretized with Chebyshev 
!                       nodes
!
!  Output arguments:
!    - npatches: integer
!        number of patches 
!    - npts: integer
!        Number of points on the discretized surface
!
!
      implicit none
      real *8, intent(in) :: rmid
      integer, intent(in) :: npars(3), norder, iptype0, iort
      integer, intent(out) :: npatches, npts

      if (iptype0.eq.1) then
        npatches = 2*(npars(1)*npars(1) + 4*npars(2)*npars(3))
        npts = npatches*(norder+1)*(norder+2)/2
      endif

      if (iptype0.eq.11.or.iptype0.eq.12) then
        npatches = npars(1)*npars(1) + 4*npars(2)*npars(3)
        npts = npatches*(norder+1)*(norder+1)
      endif

      return
      end
!
!
!
!
!
      subroutine mesh_circle_pts(rmid, npars, iort, norder, &
        iptype0, npatches, npts, ixys, ptinfo)
!
!  This subroutine returns the ptinfo on the mesh of a unit circle, 
!  in particular, it returns the u,v coordinates on points
!  on the circle and the corresponding dxy/du, and dxy/dv values  
!
!  There is middle square of length rmid, and four sectors going
!  from the vertices of the square to the boundary of the
!  circle in a radial fashion.
!
!  The middle square is discretized into (npars(1)*npars(1)) 
!  smaller squares, and each of the sectors is discretized 
!  into further (npars(2)*npars(3)) squares, with npars(2)
!  intervals in the radial direction and npars(3) 
!  intervals in the angular direction
!  
!
!  Input Arguments:
!    - rmid: real *8
!        radius of the inner square, should be between 0 and 1/sqrt(2)
!    - npars: integer(3)
!        * npars(1) = nmid, number of intervals in the middle square
!        * npars(2) = nr, number of intervals in the radial direction
!                     in each of the sectors
!        * npars(3) = nt, number of intervals in the angular direction
!                     in each of the sectors
!    - iort: integer
!        orientiation flag
!        if iort = 1, then mesh has normals pointed in +z direction
!        if iort = -1, then mesh has normals pointed in -z direction 
!    - norder: integer
!        order of discretization on each patch
!    - iptype0: integer
!        type of patch to be used in the discretization
!        * iptype = 1, triangular patch discretized using RV nodes
!        * iptype = 11, quadrangular patch discretized with GL nodes
!        * iptype = 12, quadrangular patch discretized with Chebyshev 
!                       nodes
!    - npatches: integer
!        number of patches 
!    - npts: integer
!        Number of points on the discretized surface
!
!  Output arguments:
!    - ixys: integer(npatches+1)
!        ixys(i) denotes the starting location in ptinfo
!        corresponding to patch i
!    - ptinfo: real *8 (6,npts)
!        * ptinfo(1:2,i) - x,y coordinates 
!        * ptinfo(3:4,i) - dxy/du
!        * ptinfo(5:6,i) - dxy/dv
!
!
      implicit none
      real *8, intent(in) :: rmid
      integer, intent(in) :: npars(3), norder, iptype0, iort
      integer, intent(in) :: npatches, npts

      integer, intent(out) :: ixys(npatches+1)
      real *8, intent(out) :: ptinfo(6,npts)

      real *8, allocatable :: skel(:,:,:)
      integer nskel, inds(6), ntmp, npols
      integer nmid, nr, nt
      real *8, allocatable :: uvs(:,:)
      integer i, j
      integer ipatch, ipt
      real *8 par4, dxyduv(2,2)

      procedure (), pointer :: patchpnt

      external get_xquad_circle_mesh_pt, get_xtri_circle_mesh_pt 

      nmid = npars(1)
      nr = npars(2)
      nt = npars(3)
      
      call get_npols(iptype0, norder, npols)
      allocate(uvs(2,npols))

      call get_disc_nodes(norder, npols, iptype0, uvs)
      

      allocate(skel(3,3,npatches))
      if (iptype0.eq.1) then
        call xtri_circle_mesh_skel_init(rmid, npars, iort, npatches, &
          skel, inds)
        patchpnt => get_xtri_circle_mesh_pt
      elseif (iptype0.eq.11.or.iptype0.eq.12) then
        call xquad_circle_mesh_skel_init(rmid, npars, iort, npatches, &
          skel, inds)
        patchpnt => get_xquad_circle_mesh_pt
      endif


      do ipatch = 1,npatches
        do j=1,npols
          ipt = (ipatch-1)*npols + j
          call patchpnt(ipatch, uvs(1,j), uvs(2,j), &
            ptinfo(1,ipt), dxyduv, skel, inds, rmid, par4)
          ptinfo(3,ipt) = dxyduv(1,1)
          ptinfo(4,ipt) = dxyduv(1,2)

          ptinfo(5,ipt) = dxyduv(2,1)
          ptinfo(6,ipt) = dxyduv(2,2)
        enddo
      enddo

      do i = 1,npatches+1
        ixys(i) = (i-1)*npols + 1
      enddo


      return
      end
!
!
!
!
!
      subroutine xquad_circle_mesh_skel_init(rmid, npars, iort, nskel, & 
        skel, inds)
!
!  This subroutine returns the vertices of the quadrangular skeleton mesh 
!  of a unit circle
!
!  There is middle square of length rmid, and four sectors going
!  from the vertices of the square to the boundary of the
!  circle in a radial fashion.
!
!  The middle square is discretized into (npars(1)*npars(1)) 
!  smaller squares, and each of the sectors is discretized 
!  into further (npars(2)*npars(3)) squares, with npars(2)
!  intervals in the radial direction and npars(3) 
!  intervals in the angular direction
!  
!
!  Input Arguments:
!    - rmid: real *8
!        radius of the inner square, should be between 0 and 1/sqrt(2)
!    - npars: integer(3)
!        * npars(1) = nmid, number of intervals in the middle square
!        * npars(2) = nr, number of intervals in the radial direction
!                     in each of the sectors
!        * npars(3) = nt, number of intervals in the angular direction
!                     in each of the sectors
!    - iort: integer
!        orientiation flag
!        if iort = 1, then mesh has normals pointed in +z direction
!        if iort = -1, then mesh has normals pointed in -z direction 
!   - nskel: integer
!       number of points in skeleton mesh, should be 
!       (npars(1)*npars(1) + 4*npars(2)*npars(3))
!
!  Output arguments:
!    - skel: real *8(3,3,nskel)
!        skel(1:2, ivert, ipat)
!        u,v coordinates on one of five copies of [0,1]^2, 
!        the 1st copy maps to the center region, while
!        the other four copies map to 4 sectors of width pi/2.
!
!        The third component is always 0, this is to stay
!        consistent with using xtri/xquad_rectmesh_ani
!        
!        the vertices as always are ordered as
!             3  *
!             1  2
!    - inds: integer(6)
!        inds(1): starting location of patches in the middle square
!        inds(2): starting location of patches in the first sector 
!                 (-pi/4,pi/4)
!        inds(3): starting location of patches in the second sector 
!                 (pi/4,3*pi/4)
!        inds(4): starting location of patches in the third sector 
!                 (3*pi/4,5*pi/4)
!        inds(5): starting location of patches in the fourth sector 
!                 (5*pi/4,7*pi/4)
!        inds(6): total number of skeleton patches = nskel

      implicit none
      real *8, intent(in) :: rmid
      integer, intent(in) :: npars(3), nskel
      integer, intent(in) :: iort
      real *8, intent(out) :: skel(3,3,nskel)
      integer, intent(out) :: inds(6)

      real *8 done, pi
      integer ibox, nmid, nr, nt, na, nb
      real *8 amat(2,2), vec(2)
      integer ii, nn, nout

      real *8 umin, umax, vmin, vmax
      integer nover, j
        
      done = 1
      pi   = atan(done)*4

      nmid = npars(1)
      nr = npars(2)
      nt = npars(3)

      umin = 0
      umax = 1

      vmin = 0
      vmax = 1

      if (iort.eq.-1) then
        vmin = 1
        vmax = 0
      endif

      nover = 0

      ibox = 1
      nout = 0
      inds(1) = 1
      call xquad_rectmesh_ani(umin, umax, vmin, vmax, nmid, nmid, &
        nover, nskel, nout, skel(1,1,1))
      
      ibox = ibox + nout

      do j=1,4
        inds(j+1) = ibox
        call xquad_rectmesh_ani(umin, umax, vmin, vmax, nt, nr, &
          nover, nskel, nout, skel(1,1,ibox))
        ibox = ibox + nout
      enddo
      inds(6) = ibox


      return
      end



      subroutine xtri_circle_mesh_skel_init(rmid, npars, iort, nskel, & 
        skel, inds)
!
!  This subroutine returns the vertices of the triangular skeleton mesh 
!  of a unit circle. The method proceeds by first creating a
!  quadrangular mesh and then deciding which way to split the triangles
!  based minimizing the aspect ratio
!
!  There is middle square of length rmid, and four sectors going
!  from the vertices of the square to the boundary of the
!  circle in a radial fashion.
!
!  The middle square is discretized into (npars(1)*npars(1)) 
!  smaller squares, and each of the sectors is discretized 
!  into further (npars(2)*npars(3)) squares, with npars(2)
!  intervals in the radial direction and npars(3) 
!  intervals in the angular direction
!  
!
!  Input Arguments:
!    - rmid: real *8
!        radius of the inner square, should be between 0 and 1/sqrt(2)
!    - npars: integer(3)
!        * npars(1) = nmid, number of intervals in the middle square
!        * npars(2) = nr, number of intervals in the radial direction
!                     in each of the sectors
!        * npars(3) = nt, number of intervals in the angular direction
!                     in each of the sectors
!    - iort: integer
!        orientiation flag
!        if iort = 1, then mesh has normals pointed in +z direction
!        if iort = -1, then mesh has normals pointed in -z direction 
!   - nskel: integer
!       number of points in skeleton mesh, should be 
!       (npars(1)*npars(1) + 4*npars(2)*npars(3))
!
!  Output arguments:
!    - skel: real *8(3,3,nskel)
!        skel(1:2, ivert, ipat)
!        u,v coordinates on one of five copies of [0,1]^2, 
!        the 1st copy maps to the center region, while
!        the other four copies map to 4 sectors of width pi/2.
!
!        The third component is always 0, this is to stay
!        consistent with using xtri/xquad_rectmesh_ani
!        
!        the vertices as always are ordered as
!             3  *
!             1  2
!    - inds: integer(6)
!        inds(1): starting location of patches in the middle square
!        inds(2): starting location of patches in the first sector 
!                 (-pi/4,pi/4)
!        inds(3): starting location of patches in the second sector 
!                 (pi/4,3*pi/4)
!        inds(4): starting location of patches in the third sector 
!                 (3*pi/4,5*pi/4)
!        inds(5): starting location of patches in the fourth sector 
!                 (5*pi/4,7*pi/4)
!        inds(6): total number of skeleton patches = nskel

      implicit none
      real *8, intent(in) :: rmid
      integer, intent(in) :: npars(3), nskel, iort
      real *8, intent(out) :: skel(3,3,nskel)
      integer, intent(out) :: inds(6)

      real *8 done, pi
      integer ibox, nmid, nr, nt, na, nb
      real *8 amat(2,2), vec(2)
      integer ii, nn, nout, j, nover

      real *8 umin, umax, vmin, vmax
        
      done = 1
      pi   = atan(done)*4

      nmid = npars(1)
      nr = npars(2)
      nt = npars(3)

      umin = 0
      umax = 1

      vmin = 0
      vmax = 1

      if (iort.eq.-1) then
        vmin = 1
        vmax = 0
      endif

      nover = 0

      ibox = 1
      nout = 0
      inds(1) = 1
      call xtri_rectmesh_ani(umin, umax, vmin, vmax, nmid, nmid, &
        nover, nskel, nout, skel(1,1,1))
      
      ibox = ibox + nout

      do j=1,4
        inds(j+1) = ibox
        call xquad_rectmesh_ani(umin, umax, vmin, vmax, nt, nr, &
          nover, nskel, nout, skel(1,1,ibox))
        call circ_split_quads_to_tri_minasp(nout, skel(1,1,ibox))
        ibox = ibox + 2*nout
      enddo
      inds(6) = ibox


      return
      end
!
!
!
!

      subroutine circ_split_quads_to_tri_minasp(n, skel)
      implicit real *8 (a-h,o-z)
      integer n
      real *8 skel(3,3,2*n), v1(2), v2(2), v3(2), v4(2)
      real *8, allocatable :: skeltmp(:,:,:)
      integer i, j, k

     
      
      allocate(skeltmp(3,3,n))

      do i=1,n
        do j=1,3
          do k=1,3
            skeltmp(k,j,i) = skel(k,j,i)
            skel(k,j,i) = 0
          enddo
        enddo
      enddo

      do i=1,n
         v1(1:2) = skeltmp(1:2,1,i)
         v2(1:2) = skeltmp(1:2,2,i)
         v4(1:2) = skeltmp(1:2,3,i)

         dx = v2(1) - v1(1)
         dy = v4(2) - v1(2)
         
         v3(1) = v1(1) + dx
         v3(2) = v1(2) + dy


         if ((v1(1).ge.0.5d0.and.v1(2).le.0.5d0).or. &
           (v1(1).ge.0.5d0.and.v1(2).le.0.5d0)) then
           skel(1:2,1,2*i-1) = v1(1:2)
           skel(1:2,2,2*i-1) = v2(1:2)
           skel(1:2,3,2*i-1) = v3(1:2)
           
           skel(1:2,1,2*i) = v3(1:2)
           skel(1:2,2,2*i) = v4(1:2)
           skel(1:2,3,2*i) = v1(1:2)
         else
           skel(1:2,1,2*i-1) = v4(1:2)
           skel(1:2,2,2*i-1) = v1(1:2)
           skel(1:2,3,2*i-1) = v2(1:2)
           
           skel(1:2,1,2*i) = v2(1:2)
           skel(1:2,2,2*i) = v3(1:2)
           skel(1:2,3,2*i) = v4(1:2)
         endif
      enddo

      return
      end
!
!
!
!
!
      subroutine circle_param(s, t, rmid, ang_range, xy, dxydst)
!
!  Given a point (s,t) \in [0,1]^2, this subroutine evaluates the
!  sector map, which maps the region between the line segment
!  between rmid(\cos(t0), \sin(t0)) and rmid (\cos(t1), \sin(t1))
!  and the circular arc between (\cos(t0), \sin(t0)) and
!  (\cos(t1), \sin(t1)) as a convex combination of two charts, 
!  one which maps the square to the trapezoid, and one which
!  maps the square to a cicular arc where
!  ang_range = (t0,t1)
!
      implicit real *8 (a-h,o-z)
      real *8, intent(in) :: s, t, rmid, ang_range(2)
      real *8, intent(out) :: xy(2), dxydst(2,2)
      real *8 dxy1(2,2), dxy2(2,2)
      real *8 pi
      data pi/3.141592653589793d0/

      t0 = ang_range(1)
      t1 = ang_range(2)
      
      th = (t1-t0)/2
      ta = (t1+t0)/2

      bot_wid = rmid*tan(th)
      top_wid = sin(th)

      height = cos(th) - rmid


      xt = 2*(s - 0.5d0)*(t*top_wid + (1.0d0 - t)*bot_wid)
      yt = rmid + height*t

      dxy1(1,1) = 2*(t*top_wid + (1.0d0 - t)*bot_wid)
      dxy1(1,2) = 2*(s - 0.5d0)*(top_wid - bot_wid)
      dxy1(2,1) = 0
      dxy1(2,2) = height

      thet = -th*(s-0.5d0)*2 + pi/2
      rad  = rmid/cos(th) + (1.0d0 - rmid/cos(th))*t

      xt2 = cos(thet)*rad
      yt2 = sin(thet)*rad

      dxy2(1,1) = -sin(thet)*(-th*2)*rad
      dxy2(1,2) = cos(thet)*(1.0d0 - rmid/cos(th))
      dxy2(2,1) = cos(thet)*(-th*2)*rad
      dxy2(2,2) = sin(thet)*(1.0d0 - rmid/cos(th))

      ca = cos(-ta)
      sa = sin(-ta)
      xy(1) = (xt*(1-t) + t*xt2)*ca - (yt*(1-t) + t*yt2)*sa
      xy(2) = (xt*(1-t) + t*xt2)*sa + (yt*(1-t) + t*yt2)*ca

      dxydst(1,1) = (dxy1(1,1)*(1-t) + t*dxy2(1,1))*ca - &
             (dxy1(2,1)*(1-t) + t*dxy2(2,1))*sa

      dxydst(1,2) = (dxy1(1,2)*(1-t) + t*dxy2(1,2))*ca - &
             (dxy1(2,2)*(1-t) + t*dxy2(2,2))*sa + &
              (xt2 - xt)*ca - (yt2 - yt)*sa

      dxydst(2,1) = (dxy1(1,1)*(1-t) + t*dxy2(1,1))*sa + &
             (dxy1(2,1)*(1-t) + t*dxy2(2,1))*ca
      
      dxydst(2,2) = (dxy1(1,2)*(1-t) + t*dxy2(1,2))*sa + &
             (dxy1(2,2)*(1-t) + t*dxy2(2,2))*ca + &
             (xt2- xt)*sa + (yt2-yt)*ca
      

      return
      end
!
!
!
!
!
!
      subroutine get_xtri_circle_mesh_pt(itri, u, v, xy, dxyduv, &
        skel, ipars, rmid, par4)
      
!
!  Given a triangle id on the skeleton mesh, this routine 
!  first identifies which of the 5 faces of the skeleton mesh
!  the triangle is on, and then appropriately discretizes the
!  point on the circle mesh.
!
!  ipars(1:6) = inds(1:6)
!
      implicit none
      integer, intent(in) :: itri
      real *8, intent(in) :: u, v, skel(3,3,*)
      integer, intent(in) :: ipars(6)
      real *8, intent(in) :: rmid
      real *8 par4
      real *8, intent(out) :: xy(2), dxyduv(2,2)

      integer iface
      real *8 x0, y0, x1, y1, x2, y2
      real *8 dsdu, dtdu, dsdv, dtdv, s, t, dxydst(2,2)
      real *8 ang_range(2)
      real *8 pi
      data pi/3.141592653589793d0/

      

      if (itri.lt.ipars(2)) then
        iface = 0
      elseif (itri.lt.ipars(3)) then
        iface = 1
      elseif (itri.lt.ipars(4)) then
        iface = 2
      elseif (itri.lt.ipars(5)) then
        iface = 3
      elseif (itri.lt.ipars(6)) then
        iface = 4
      endif

      x0 = skel(1,1,itri)
      y0 = skel(2,1,itri)

      x1 = skel(1,2,itri)
      y1 = skel(2,2,itri)

      x2 = skel(1,3,itri)
      y2 = skel(2,3,itri)
        
      dsdu = (x1-x0)
      dsdv = (x2-x0)

      dtdu = (y1-y0)
      dtdv = (y2-y0)

      s = x0 + (x1-x0)*u + (x2-x0)*v
      t = y0 + (y1-y0)*u + (y2-y0)*v

      ang_range(1) = -pi/4 + pi/2*(iface-1)
      ang_range(2) =  pi/4 + pi/2*(iface-1)
      if (iface.eq.0) then
        xy(1) = 2*rmid*s - rmid 
        xy(2) = 2*rmid*t - rmid

        dxyduv(1,1) = 2*rmid*dsdu
        dxyduv(1,2) = 2*rmid*dtdu

        dxyduv(2,1) = 2*rmid*dsdv
        dxyduv(2,2) = 2*rmid*dtdv
      else
        call circle_param(s, t, rmid, ang_range, xy, dxydst)
            
        dxyduv(1,1) = dxydst(1,1)*dsdu + dxydst(1,2)*dtdu
        dxyduv(1,2) = dxydst(2,1)*dsdu + dxydst(2,2)*dtdu

        dxyduv(2,1) = dxydst(1,1)*dsdv + dxydst(1,2)*dtdv
        dxyduv(2,2) = dxydst(2,1)*dsdv + dxydst(2,2)*dtdv
      endif

      return
      end
!
!
!
!
      subroutine get_xquad_circle_mesh_pt(iquad, u, v, xy, dxyduv, &
        skel, ipars, rmid, par4)
      
!
!  Given a quad id on the skeleton mesh, this routine 
!  first identifies which of the 5 faces of the skeleton mesh
!  the triangle is on, and then appropriately discretizes the
!  point on the circle mesh.
!
!  ipars(1:6) = inds(1:6)
!
      implicit none
      integer, intent(in) :: iquad
      real *8, intent(in) :: u, v, skel(3,3,*)
      integer, intent(in) :: ipars(6)
      real *8, intent(in) :: rmid
      real *8 par4
      real *8, intent(out) :: xy(2), dxyduv(2,2)

      integer iface
      real *8 x0, y0, x1, y1, x2, y2
      real *8 dsdu, dtdu, dsdv, dtdv, s, t, dxydst(2,2)
      real *8 ang_range(2)
      real *8 pi
      data pi/3.141592653589793d0/

      

      if (iquad.lt.ipars(2)) then
        iface = 0
      elseif (iquad.lt.ipars(3)) then
        iface = 1
      elseif (iquad.lt.ipars(4)) then
        iface = 2
      elseif (iquad.lt.ipars(5)) then
        iface = 3
      elseif (iquad.lt.ipars(6)) then
        iface = 4
      endif

      x0 = skel(1,1,iquad)
      y0 = skel(2,1,iquad)

      x1 = skel(1,2,iquad)
      y1 = skel(2,2,iquad)

      x2 = skel(1,3,iquad)
      y2 = skel(2,3,iquad)
        
      dsdu = (x1-x0)/2
      dsdv = (x2-x0)/2

      dtdu = (y1-y0)/2
      dtdv = (y2-y0)/2

      s = x0 + (x1-x0)*(1.0d0+u)/2 + (x2-x0)*(1.0d0+v)/2
      t = y0 + (y1-y0)*(1.0d0+u)/2 + (y2-y0)*(1.0d0+v)/2

      ang_range(1) = -pi/4 + pi/2*(iface-1)
      ang_range(2) =  pi/4 + pi/2*(iface-1)
      if (iface.eq.0) then
        xy(1) = 2*rmid*s - rmid 
        xy(2) = 2*rmid*t - rmid

        dxyduv(1,1) = 2*rmid*dsdu
        dxyduv(1,2) = 2*rmid*dtdu

        dxyduv(2,1) = 2*rmid*dsdv
        dxyduv(2,2) = 2*rmid*dtdv
      else
        call circle_param(s, t, rmid, ang_range, xy, dxydst)
            
        dxyduv(1,1) = dxydst(1,1)*dsdu + dxydst(1,2)*dtdu
        dxyduv(1,2) = dxydst(2,1)*dsdu + dxydst(2,2)*dtdu

        dxyduv(2,1) = dxydst(1,1)*dsdv + dxydst(1,2)*dtdv
        dxyduv(2,2) = dxydst(2,1)*dsdv + dxydst(2,2)*dtdv
      endif

      return
      end
