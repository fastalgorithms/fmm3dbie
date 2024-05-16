c     This routines assumes the closest point is a surface discretization point
c     flag = 2: inside
c     flag = 1: outside
   
      subroutine inoutsurf3d(npatches,norder,npols,nptssurf,
     1     srccoefs,ntargs,targs,sxyz,itvec,uvsvec,flags)
      
      implicit none
      
c     List of calling arguments
      integer npatches,norder,nptssurf,npols,flags(ntargs)
      integer itvec(ntargs),ntargs
      real *8 srccoefs(9,nptssurf)
      real *8 uvsvec(2,ntargs),targs(3,ntargs)
      
c     List of local variables
      integer ipnt,k,npols2,it
      real *8 rdotn,sxyz(3,ntargs),nxyz(3),r(3)
      real *8 dxyzdu(3),dxyzdv(3),pols(2000)
            
      do ipnt = 1,ntargs
c     starting index for the patch with the closest surface
c     point to targs(:,ipnt)
         it = itvec(ipnt)
       
c     Convert (u,v) for the closest point, in local parameter 
c     space, to a point in R3^
        
         call koorn_pols(uvsvec(1,ipnt),norder,npols2,pols)
         dxyzdu(1) = 0.0d0
         dxyzdu(2) = 0.0d0
         dxyzdu(3) = 0.0d0
         dxyzdv(1) = 0.0d0
         dxyzdv(2) = 0.0d0
         dxyzdv(3) = 0.0d0

c     interpolate derivatives wrt local parametrization
         
         do k=1,npols
            dxyzdu(1) = dxyzdu(1)
     1           + srccoefs(4,it+k-1)*pols(k)
            dxyzdu(2) = dxyzdu(2)
     1           + srccoefs(5,it+k-1)*pols(k)
            dxyzdu(3) = dxyzdu(3)
     1           + srccoefs(6,it+k-1)*pols(k)
            dxyzdv(1) = dxyzdv(1)
     1           + srccoefs(7,it+k-1)*pols(k)
            dxyzdv(2) = dxyzdv(2)
     1           + srccoefs(8,it+k-1)*pols(k)
            dxyzdv(3) = dxyzdv(3)
     1           + srccoefs(9,it+k-1)*pols(k)
         end do
         
         r(1) = targs(1,ipnt)-sxyz(1,ipnt)
         r(2) = targs(2,ipnt)-sxyz(2,ipnt)
         r(3) = targs(3,ipnt)-sxyz(3,ipnt)
c     get normal at closest surface point
         call cross_prod3d(dxyzdu,dxyzdv,nxyz)
c     inner product with r and normal
         rdotn = r(1)*nxyz(1)+r(2)*nxyz(2)+r(3)*nxyz(3)
c     check normal sign
         flags(ipnt) = -1
c     inside
         if(rdotn.lt.0.0d0) then
            flags(ipnt) = 2
         endif
c     outside
         if(rdotn.gt.0.0d0) then
            flags(ipnt) = 1
         endif
c         
      enddo      

      return
      end
c     end of subroutine inoutsurf3d
c     
c
c
