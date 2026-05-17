c     
c     given a target point, this routine returns the closest point on a patch, using Newtons method,
c     the initial guess must be a surface discretization point.
c     
      
      subroutine findnearsrfcpnt(ntargs,targs,
     1     nptssurf,porder,snpols,maxiter,srcvals,srccoefs,tol,
     2     lpatchidxvec,itvec,sxyz,uvsloc,dists,flags)
      
      implicit none

c     List of calling arguments
      integer ntargs,flags(ntargs)
      integer nptssurf,porder,maxiter
      integer itvec(ntargs),lpatchidxvec(ntargs)
      
      real *8 srcvals(12,nptssurf),srccoefs(9,nptssurf)
      real *8 targs(3,ntargs),sxyz(3,ntargs)
      real *8 dists(ntargs),tol,uvsloc(2,ntargs)

c     List of local variables
      integer snpols,snpols2,nn,it,lpatchidx
      integer ipnt,idx,k
      real *8 sxyz0(3),bs,uv0(2),pols(2000)
      real *8 uv_nodes(2,snpols),umatr(snpols,snpols)
      real *8 vmatr(snpols,snpols),wts(snpols)
      real *8 rsc1(0:porder, 0:porder), rat1(2, 0:porder)
      real *8 rat2(3, 0:porder, 0:porder)

c     get Vioreanu Rokhlin nodes on simplex
      call vioreanu_simplex_quad(porder,snpols,uv_nodes,umatr,vmatr,wts)
c     precompute Koornwinder polys
      call koornf_init(porder, rat1, rat2, rsc1)
c     loop over target points
      do ipnt = 1,ntargs
         lpatchidx = lpatchidxvec(ipnt)
         it = itvec(ipnt)
c     Get an initial guess in (u,v) parameter space      
         uv0(1) = uv_nodes(1,lpatchidx)
         uv0(2) = uv_nodes(2,lpatchidx)
c     Newton to find closest point, returns point (u,v) in parameter space
         call newton_nearp_fast(srcvals,srccoefs,nptssurf,it,
     1        porder,targs(1,ipnt),rat1,rat2,rsc1,uv0,uvsloc(1,ipnt),
     2        tol,maxiter,flags(ipnt))
c     Convert (u,v) in parameter space to a point sxyz in R^3
         call koorn_pols(uvsloc(1,ipnt),porder,snpols2, pols)
         sxyz(1,ipnt) = 0.0d0
         sxyz(2,ipnt) = 0.0d0
         sxyz(3,ipnt) = 0.0d0        
         do k=1,snpols
            sxyz(1,ipnt) = sxyz(1,ipnt) + srccoefs(1,it+k-1)*pols(k)
            sxyz(2,ipnt) = sxyz(2,ipnt) + srccoefs(2,it+k-1)*pols(k)
            sxyz(3,ipnt) = sxyz(3,ipnt) + srccoefs(3,it+k-1)*pols(k)
         end do
         
c     Compute distance between target point its closest point sxyz
         dists(ipnt)=(sxyz(1,ipnt)-targs(1,ipnt))**2
     1        + (sxyz(2,ipnt)-targs(2,ipnt))**2
     2        + (sxyz(3,ipnt)-targs(3,ipnt))**2
         dists(ipnt) = dsqrt(dists(ipnt))
      enddo
      return    
c     end of findnearsrfcpnt
      end 
c     
c     
      


      subroutine newton_nearp(srcvals,srccoefs,npoints,it,order,pt,uv0,
     1     uvs,tol,maxiter,flag)
      implicit none

!List of calling arguments
      integer ( kind = 4 ), intent(in) :: npoints
      integer ( kind = 4 ), intent(in) :: order
      integer ( kind = 4 ), intent(in) :: it,maxiter
      real ( kind = 8 ), intent(in) :: srcvals(12,npoints)
      real ( kind = 8 ), intent(in) :: srccoefs(9,npoints)
      real ( kind = 8 ), intent(in) :: pt(3)
      real ( kind = 8 ), intent(in) :: tol
      real ( kind = 8 ), intent(inout) :: uv0(2)
      real ( kind = 8 ), intent(out) :: uvs(2)
      integer ( kind = 4 ), intent(out) :: flag
      

!List of local variables
      real ( kind = 8 ) err,aux2
      integer ( kind = 4 ) iter
      real ( kind = 8 ) gradf(2),hesf(2,2),delta(2),dethesf

      err=(tol+1.0d0)**2
      iter=0

      do while ((iter<maxiter).and.(err>tol**2))

         call eval_fgradhess(srcvals,srccoefs,npoints,it,order,uv0,pt,
     1        gradf,hesf)
         dethesf=hesf(1,1)*hesf(2,2)-hesf(1,2)*hesf(2,1)
         delta(1)=-(hesf(1,2)*gradf(2)-hesf(2,2)*gradf(1))/dethesf
         delta(2)=(hesf(1,1)*gradf(2)-hesf(2,1)*gradf(1))/dethesf
         

         aux2=delta(1)*gradf(1)+delta(2)*gradf(2) 
         if (aux2<0.0d0) then
            delta(1)=0.1d0*gradf(1)
            delta(2)=0.1d0*gradf(2)          
         endif

         uvs(1)=uv0(1)-delta(1)
         uvs(2)=uv0(2)-delta(2)

         iter=iter+1

         if ((uvs(1)+uvs(2))>1.0d0) then
            uvs=uvs/(uvs(1)+uvs(2))
         endif
         
         if (uvs(1)<0.0d0) then
            uvs(1)=0.0d0
         endif

         if (uvs(2)<0.0d0) then
            uvs(2)=0.0d0
         endif

         if (uvs(1)>1.0d0) then
            uvs(1)=1.0d0
         endif

         if (uvs(2)>1.0d0) then
            uvs(2)=1.0d0
         endif
         
         err=(uvs(1)-uv0(1))**2+(uvs(2)-uv0(2))**2
         
         uv0(1)=uvs(1)
         uv0(2)=uvs(2)

      enddo

      if (err<tol**2) then
         flag=0
      else
         flag=1
      endif

      return
      end subroutine newton_nearp
!     
!     
!     
!     
!     

      subroutine newton_nearp_fast(srcvals,srccoefs,npoints,it,order,pt,
     1     rat1, rat2, rsc1, uv0,uvs,tol,maxiter,flag)
      implicit none

!List of calling arguments
      integer ( kind = 4 ), intent(in) :: npoints
      integer ( kind = 4 ), intent(in) :: order
      integer ( kind = 4 ), intent(in) :: it,maxiter
      real ( kind = 8 ), intent(in) :: srcvals(12,npoints)
      real ( kind = 8 ), intent(in) :: srccoefs(9,npoints)
      real ( kind = 8 ), intent(in) :: pt(3)
      real ( kind = 8 ), intent(in) :: tol,rat1(2,0:order)
      real ( kind = 8 ), intent(in) :: rsc1(0:order,0:order)
      real ( kind = 8 ), intent(in) :: rat2(3, 0:order, 0:order)
      real ( kind = 8 ), intent(inout) :: uv0(2)
      real ( kind = 8 ), intent(out) :: uvs(2)
      integer ( kind = 4 ), intent(out) :: flag
      

!List of local variables
      real ( kind = 8 ) err,aux2,tritol,tritolu,tritoll
      integer ( kind = 4 ) iter
      real ( kind = 8 ) gradf(2),hesf(2,2),delta(2),dethesf

      err=(tol+1.0d0)**2
      iter=0

      do while ((iter<maxiter).and.(err>tol**2))

         call eval_fgradhess_fast(srcvals,srccoefs,npoints,it,order,uv0, 
     1        rat1, rat2, rsc1, pt,gradf,hesf)
         dethesf=hesf(1,1)*hesf(2,2)-hesf(1,2)*hesf(2,1)
         delta(1)=-(hesf(1,2)*gradf(2)-hesf(2,2)*gradf(1))/dethesf
         delta(2)=(hesf(1,1)*gradf(2)-hesf(2,1)*gradf(1))/dethesf
         

         aux2=delta(1)*gradf(1)+delta(2)*gradf(2) 
         if (aux2<0.0d0) then
            delta(1)=0.1d0*gradf(1)
            delta(2)=0.1d0*gradf(2)          
         endif

         uvs(1)=uv0(1)-delta(1)
         uvs(2)=uv0(2)-delta(2)

         iter=iter+1

         tritol = 1.0d-6
         tritolu = 1.0d0+tritol
         tritoll = 0.0d0-tritol
c     if ((uvs(1)+uvs(2))>1.0d0) then
         if ((uvs(1)+uvs(2))>tritolu) then
            uvs=uvs/(uvs(1)+uvs(2))
         endif
         
c         if (uvs(1)<0.0d0) then
         if (uvs(1)<tritoll) then
            uvs(1)=tritoll
         endif

         if (uvs(2)<tritoll) then
            uvs(2)=tritoll
         endif

         if (uvs(1)>tritolu) then
            uvs(1)=tritolu
         endif

         if (uvs(2)>tritolu) then
            uvs(2)=tritolu
         endif
         
         err=(uvs(1)-uv0(1))**2+(uvs(2)-uv0(2))**2
         
         uv0(1)=uvs(1)
         uv0(2)=uvs(2)

      enddo

      if (err<tol**2) then
         flag=0
      else
         flag=1
      endif

      return
      end subroutine newton_nearp_fast


!     
!     
!     
!     
      subroutine eval_fgradhess(srcvals,srccoefs,npoints,it,order,uvs,
     1     pt,gradf,hesf)
      implicit none

!List of calling arguments
      integer ( kind = 4 ), intent(in) :: npoints
      integer ( kind = 4 ), intent(in) :: order
      integer ( kind = 4 ), intent(in) :: it
      real ( kind = 8 ), intent(in) :: srcvals(12,npoints)
      real ( kind = 8 ), intent(in) :: srccoefs(9,npoints)
      real ( kind = 8 ), intent(in) :: pt(3)
      real ( kind = 8 ), intent(in) :: uvs(2)
      real ( kind =  8 ), intent(out) :: gradf(2)
      real ( kind = 8 ), intent(out) :: hesf(2,2)
      
!List of local variables
      integer ( kind = 4 ) npols, i 
      real ( kind = 8 ) x,y,z,xu,yu,zu,xv,yv,zv
      real ( kind = 8 ) xuu,xuv,xvv,yuu,yuv,yvv,zuu,zuv,zvv,aux1(2)

!     Note that these could be made allocatable if needed
      real *8 pols(500), ders(2,500)
      

      npols = (order+1)*(order+2)/2

      call koorn_ders(uvs, order, npols, pols, ders)
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
         x = x + srccoefs(1,it+i-1)*pols(i)
         y = y + srccoefs(2,it+i-1)*pols(i)
         z = z + srccoefs(3,it+i-1)*pols(i)
         
         xu = xu + srccoefs(4,it+i-1)*pols(i)
         yu = yu + srccoefs(5,it+i-1)*pols(i)
         zu = zu + srccoefs(6,it+i-1)*pols(i)

         xv = xv + srccoefs(7,it+i-1)*pols(i)
         yv = yv + srccoefs(8,it+i-1)*pols(i)
         zv = zv + srccoefs(9,it+i-1)*pols(i)

         xuu = xuu + srccoefs(4,it+i-1)*ders(1,i)
         yuu = yuu + srccoefs(5,it+i-1)*ders(1,i)
         zuu = zuu + srccoefs(6,it+i-1)*ders(1,i)

         xuv = xuv + srccoefs(4,it+i-1)*ders(2,i)
         yuv = yuv + srccoefs(5,it+i-1)*ders(2,i)
         zuv = zuv + srccoefs(6,it+i-1)*ders(2,i)

         xvv = xvv + srccoefs(7,it+i-1)*ders(2,i)
         yvv = yvv + srccoefs(8,it+i-1)*ders(2,i)
         zvv = zvv + srccoefs(9,it+i-1)*ders(2,i)
      enddo


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
      
      subroutine eval_fgradhess_fast(srcvals,srccoefs,npoints,it,
     1     order,uvs,rat1,rat2,rsc1, pt,gradf,hesf)
!     
!     Given srcvals, and srccoefs, evalute, xyz, gradient and hessian info
!     at a new point. This is the fast version which includes 
!     precomputed arrays for evaluating the gradient and hessian info
!     

      implicit none

!List of calling arguments
      integer ( kind = 4 ), intent(in) :: npoints
      integer ( kind = 4 ), intent(in) :: order
      integer ( kind = 4 ), intent(in) :: it
      real ( kind = 8 ), intent(in) :: srcvals(12,npoints)
      real ( kind = 8 ), intent(in) :: srccoefs(9,npoints)
      real ( kind = 8 ), intent(in) :: pt(3)
      real ( kind = 8 ), intent(in) :: uvs(2),rat1(2,0:order)
      real ( kind = 8 ), intent(in) :: rsc1(0:order,0:order)
      real ( kind = 8 ), intent(in) :: rat2(3, 0:order, 0:order)
      real ( kind =  8 ), intent(out) :: gradf(2)
      real ( kind = 8 ), intent(out) :: hesf(2,2)
      
!List of local variables
      integer ( kind = 4 ) npols, i
      real ( kind = 8 ) x,y,z,xu,yu,zu,xv,yv,zv
      real ( kind = 8 ) xuu,xuv,xvv,yuu,yuv,yvv,zuu,zuv,zvv,aux1(2)

!     Note that these could be made allocatable if needed
      real *8 pols(500), ders(2,500)
      

      npols = (order+1)*(order+2)/2
      call koornf_ders(uvs,order, npols, pols, ders, rat1, rat2, rsc1)
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
         x = x + srccoefs(1,it+i-1)*pols(i)
         y = y + srccoefs(2,it+i-1)*pols(i)
         z = z + srccoefs(3,it+i-1)*pols(i)
         
         xu = xu + srccoefs(4,it+i-1)*pols(i)
         yu = yu + srccoefs(5,it+i-1)*pols(i)
         zu = zu + srccoefs(6,it+i-1)*pols(i)

         xv = xv + srccoefs(7,it+i-1)*pols(i)
         yv = yv + srccoefs(8,it+i-1)*pols(i)
         zv = zv + srccoefs(9,it+i-1)*pols(i)

         xuu = xuu + srccoefs(4,it+i-1)*ders(1,i)
         yuu = yuu + srccoefs(5,it+i-1)*ders(1,i)
         zuu = zuu + srccoefs(6,it+i-1)*ders(1,i)

         xuv = xuv + srccoefs(4,it+i-1)*ders(2,i)
         yuv = yuv + srccoefs(5,it+i-1)*ders(2,i)
         zuv = zuv + srccoefs(6,it+i-1)*ders(2,i)

         xvv = xvv + srccoefs(7,it+i-1)*ders(2,i)
         yvv = yvv + srccoefs(8,it+i-1)*ders(2,i)
         zvv = zvv + srccoefs(9,it+i-1)*ders(2,i)
      enddo


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
      end subroutine eval_fgradhess_fast




