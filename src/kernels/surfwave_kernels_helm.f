! Kernels for Helmholtz surface wave problems
!
! the following routines rely on the srcinfo and targinfo arrays
! containing the following information, standardized in the following
! way:
!
!     *info(1:3) = xyz
!     *info(4:5) = tangent vector in the 2d plane
!     *info(10:11) = normal vector in the 2d plane
!
!     *zpars(1:3) = roots of the cubic dispersion polynomial,
!           zpars(1) must be the root closest to the real axis
!     *zpars(4:6) = residues corresponding to zpars(1:3)
!
! The Helmholtz surface wave Green's functions gs and gphi are sums
! over the three dispersion roots, evaluated via the primitive routines
! gshelm and gphihelm. Gradient routines operate on (dx,dy) directly.
!
c
c
c  Kernel wrappers (src/targ in srcinfo/targinfo format)
c
c
        subroutine gshelmkern(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,
     1    ipars,val)
c
c  returns the Helmholtz surface wave free-surface Green's function
c  G_s(r) at r = |x-y|
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz,dr
        complex *16 zpars(ndz), val        
        complex *16 rts(3), ejs(3)

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)

        dr=sqrt(dx**2+dy**2+dz**2)

        rts = zpars(1:3)
        ejs = zpars(4:6) 

        call gshelm(rts,ejs,dr,val)

        return 
        end
c
c
c
        subroutine gshelm_sp_kern(src,ndt,targ,ndd,dpars,ndz,zpars,
     1       ndi,ipars,val)
c
c  returns the normal derivative of G_s at x-y = (dx,dy),
c  dotted with the target normal targ(10:11):
c    val = grad G_s . n
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz,dr, nx, ny
        complex *16 zpars(ndz), val        
        complex *16 rts(3), ejs(3), gsx, gsy

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)

        dr=sqrt(dx**2+dy**2+dz**2)

        rts = zpars(1:3)
        ejs = zpars(4:6)

        nx = targ(10)
        ny = targ(11)

        call gradxgshelm(rts,ejs,dx,dy,gsx)
        call gradygshelm(rts,ejs,dx,dy,gsy)

        val = gsx*nx+gsy*ny

        return 
        end
c
c
c
c
        subroutine gphihelm_sp_kern(src,ndt,targ,ndd,dpars,ndz,zpars,
     1       ndi,ipars,val)
c
c  returns the normal derivative of the velocity potential Green's
c  function G_phi, dotted with the target normal targ(10:11):
c    val = grad G_phi . n
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz,dr, nx, ny
        complex *16 zpars(ndz), val
        complex *16 rts(3), ejs(3), gsx, gsy

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)

        dr=sqrt(dx**2+dy**2+dz**2)

        rts = zpars(1:3)
        ejs = zpars(4:6)

        nx = targ(10)
        ny = targ(11)

        call gradxgphihelm(rts,ejs,dx,dy,gsx)
        call gradygphihelm(rts,ejs,dx,dy,gsy)

        val = gsx*nx+gsy*ny

        return
        end
c
c
c
c
        subroutine gphihelmkern(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,
     1    ipars,val)
c
c  returns the Helmholtz surface wave velocity potential Green's function
c  G_phi(r) at r = |x-y|
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz,dr
        complex *16 zpars(ndz), val
        complex *16 rts(3), ejs(3)

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)

        dr=sqrt(dx**2+dy**2+dz**2)

        rts = zpars(1:3)
        ejs = zpars(4:6)

        call gphihelm(rts,ejs,dr,val)

        return
        end

c
c
c
c
        subroutine s3dgphihelmkern(src,ndt,targ,ndd,dpars,ndz,zpars,
     1    ndi,ipars,val)
c
c  returns the 3D free-space contribution to G_phi at r = |x-y|
c  (residue-weighted sum without the rts prefactor, scaled by 1/8)
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz,dr
        complex *16 zpars(ndz), val
        complex *16 rts(3), ejs(3)

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)

        dr=sqrt(dx**2+dy**2+dz**2)

        rts = zpars(1:3)
        ejs = zpars(4:6)

        call s3dgphihelm(rts,ejs,dr,val)

        return
        end

c
c
c
c
        subroutine lapgshelmkern(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,
     1    ipars,val)
c
c  returns the Laplacian of the free-surface Green's function: Delta G_s(r)
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz,dr
        complex *16 zpars(ndz), val
        complex *16 rts(3), ejs(3)

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)

        dr=sqrt(dx**2+dy**2+dz**2)

        rts = zpars(1:3)
        ejs = zpars(4:6)

        call lapgshelm(rts,ejs,dr,val)

        return
        end
c
c
c
c
c
        subroutine lapgphihelmkern(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,
     1    ipars,val)
c
c  returns the Laplacian of the velocity potential Green's function:
c  Delta G_phi(r), including the 1/r free-space correction
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz,dr
        complex *16 zpars(ndz), val
        complex *16 rts(3), ejs(3)

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)

        dr=sqrt(dx**2+dy**2+dz**2)

        rts = zpars(1:3)
        ejs = zpars(4:6)

        call lapgphihelm(rts,ejs,dr,val)

        return
        end
c
c
c
c
        subroutine gradxgshelmkern(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,
     1    ipars,val)
c
c  returns the x-derivative of the free-surface Green's function:
c  d G_s / dx at (dx,dy) = targ(1:2) - src(1:2)
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz,dr
        complex *16 zpars(ndz), val
        complex *16 rts(3), ejs(3)

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)

        dr=sqrt(dx**2+dy**2+dz**2)

        rts = zpars(1:3)
        ejs = zpars(4:6)

        call gradxgshelm(rts,ejs,dx,dy,val)

        return
        end
c
c
c
c
        subroutine gradygshelmkern(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,
     1    ipars,val)
c
c  returns the y-derivative of the free-surface Green's function:
c  d G_s / dy at (dx,dy) = targ(1:2) - src(1:2)
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz,dr
        complex *16 zpars(ndz), val
        complex *16 rts(3), ejs(3)

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)      
        dz=targ(3)-src(3)

        dr=sqrt(dx**2+dy**2+dz**2)

        rts = zpars(1:3)
        ejs = zpars(4:6)

        call gradygshelm(rts,ejs,dx,dy,val)

        return
        end
c
c
c
        subroutine gradxgphihelmkern(src,ndt,targ,ndd,dpars,ndz,zpars,
     1    ndi,ipars,val)
c
c  returns the x-derivative of the velocity potential Green's function:
c  d G_phi / dx at (dx,dy) = targ(1:2) - src(1:2)
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz,dr
        complex *16 zpars(ndz), val
        complex *16 rts(3), ejs(3)

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)

        dr=sqrt(dx**2+dy**2+dz**2)

        rts = zpars(1:3)
        ejs = zpars(4:6)

        call gradxgphihelm(rts,ejs,dx,dy,val)

        return
        end
c
c
c
c
        subroutine gradygphihelmkern(src,ndt,targ,ndd,dpars,ndz,zpars,
     1    ndi,ipars,val)
c
c  returns the y-derivative of the velocity potential Green's function:
c  d G_phi / dy at (dx,dy) = targ(1:2) - src(1:2)
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz,dr
        complex *16 zpars(ndz), val
        complex *16 rts(3), ejs(3)

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)

        dr=sqrt(dx**2+dy**2+dz**2)

        rts = zpars(1:3)
        ejs = zpars(4:6)

        call gradygphihelm(rts,ejs,dx,dy,val)

        return
        end
c
c

c
c
c
c
c
c
c
c  Primitive Green's function evaluators (r-only or (dx,dy) interface)
c
        subroutine gshelm(rts,ejs,dr,val)
c
c  returns the free-surface Green's function G_s(r) as a weighted sum
c  over 3 dispersion roots: real root via Hankel+K_0, complex roots via K_0
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(3),ejs(3),val,zt
        complex *16 h0,h1,sk01,sk02,sk03,ima,rhoj
        real *8 dr,pi
        logical ilow
        data ima /(0,1)/

        rhoj = rts(1)
       
        if (abs(atan(imag(rhoj) / real(rhoj))) .le. 1d-14) then 
             zt = rts(1)*dr
        ione = 1
             call hank103(zt,h0,h1,ione)
             call sk0(rhoj,dr,sk01)
             val = ejs(1)*rts(1)**2*(-sk01+2*ima*h0)
        endif
        if (abs(atan(imag(rhoj) / real(rhoj))) .ge. 1d-14) then
             rhoj = -rhoj
             call sk0(rhoj,dr,sk01)
             val = ejs(1)*rts(1)**2*sk01
        endif
        
        rhoj = -rts(2)
        call sk0(rhoj,dr,sk02)

        val = val + ejs(2)*rts(2)**2*sk02
        
        rhoj = -rts(3)
        call sk0(rhoj,dr,sk03)

        val = val + ejs(3)*rts(3)**2*sk03
        val = val/2

        return
        end
c
c
c
        subroutine lapgshelm(rts,ejs,dr,val)
c
c  returns the Laplacian of G_s: Delta G_s(r) = -sum_j ejs(j)*rts(j)^4 * K_0(-rts(j)*r)
c  (real root uses Hankel, complex roots use K_0)
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(3),ejs(3),val,zt
        complex *16 h0,h1,sk01,sk02,sk03,ima,rhoj
        real *8 dr,pi
        logical ilow
        data ima /(0,1)/

        zt = rts(1)*dr

        rhoj = rts(1)

        ione = 1
        call hank103(zt,h0,h1,ione)
        call sk0(rhoj,dr,sk01)

        val = ejs(1)*rts(1)**4*(sk01-2*ima*h0)

        rhoj = -rts(2)
        call sk0(rhoj,dr,sk02)

        val = val - ejs(2)*rts(2)**4*sk02

        rhoj = -rts(3)
        call sk0(rhoj,dr,sk03)

        val = val - ejs(3)*rts(3)**4*sk03
        val = val/2

        return
        end
c
c
c
c
c
        subroutine gphihelm(rts,ejs,dr,val)
c
c  returns the velocity potential Green's function G_phi(r): same structure
c  as gshelm but with rts(j)^1 prefactor instead of rts(j)^2, scaled by 1/4
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(3),ejs(3),val,zt
        complex *16 h0,h1,sk01,sk02,sk03,ima,rhoj
        real *8 dr,pi
        logical ilow
        data ima /(0,1)/

        rhoj = rts(1)

        if (abs(atan(imag(rhoj) / real(rhoj))) .le. 1d-14) then
             zt = rts(1)*dr
        ione = 1
             call hank103(zt,h0,h1,ione)
             call sk0(rhoj,dr,sk01)
             val = ejs(1)*rts(1)*(-sk01+2*ima*h0)
        endif
        if (abs(atan(imag(rhoj) / real(rhoj))) .ge. 1d-14) then
             rhoj = -rhoj
             call sk0(rhoj,dr,sk01)
             val = ejs(1)*rts(1)*sk01
        endif

        rhoj = -rts(2)
        call sk0(rhoj,dr,sk02)

        val = val + ejs(2)*rts(2)*sk02

        rhoj = -rts(3)
        call sk0(rhoj,dr,sk03)

        val = val + ejs(3)*rts(3)*sk03
        val = val/4

        return
        end
c
c
c
        subroutine s3dgphihelm(rts,ejs,dr,val)
c
c  returns the 3D free-space part of G_phi: residue-weighted K_0 sum
c  without the rts prefactor, scaled by 1/8
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(3),ejs(3),val,zt
        complex *16 h0,h1,sk01,sk02,sk03,ima,rhoj
        real *8 dr,pi
        logical ilow
        data ima /(0,1)/
        
        rhoj = rts(1)
        
        if (abs(atan(imag(rhoj) / real(rhoj))) .le. 1d-14) then
             zt = rts(1)*dr
        ione = 1
             call hank103(zt,h0,h1,ione)
             call sk0(rhoj,dr,sk01)
             val = ejs(1)*(-sk01+2*ima*h0)
        endif
        if (abs(atan(imag(rhoj) / real(rhoj))) .ge. 1d-14) then
             rhoj = -rhoj
             call sk0(rhoj,dr,sk01)
             val = ejs(1)*sk01
        endif

        rhoj = -rts(2)
        call sk0(rhoj,dr,sk02)
        
        val = val + ejs(2)*sk02
        
        rhoj = -rts(3)
        call sk0(rhoj,dr,sk03)

        val = val + ejs(3)*sk03
        val = val/8

        return
        end
c
c
c
c
        subroutine lapgphihelm_old(rts,ejs,dr,val)
c
c  (deprecated) earlier version of lapgphihelm using lapsk0 directly
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(3),ejs(3),val,zt,laph0
        complex *16 h0,h1,lap1,lap2,lap3,ima,rhoj
        real *8 dr,pi
        logical ilow  
        data ima /(0,1)/

        zt = rts(1)*dr

        rhoj = rts(1)

        ione = 1
        call hank103(zt,h0,h1,ione)
        call lapsk0(rhoj,dr,lap1)

        laph0 = -rhoj*rhoj*h0
        val = ejs(1)*rts(1)*(-lap1+2*ima*laph0)

        rhoj = -rts(2)
        call lapsk0(rhoj,dr,lap2)

        val = val + ejs(2)*rts(2)*(lap2)

        rhoj = -rts(3)   
        call lapsk0(rhoj,dr,lap3)

        val = val + ejs(3)*rts(3)*(lap3)
        val = val/4

        return
        end
c
c
c
c
        subroutine lapgphihelm(rts,ejs,dr,val)
c
c  returns Delta G_phi(r), including the 1/r free-space correction term
c  (uses alpha = 1/sum_j ejs(j)*rts(j)^2 to scale the correction)
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(3),ejs(3),val,zt,alpha
        complex *16 h0,h1,sk01,sk02,sk03,ima,rhoj
        real *8 dr,pi
        logical ilow
        data ima /(0,1)/,pi/0.31415926535897932D+01/

        zt = rts(1)*dr
        rhoj = rts(1)

        ione = 1
        call hank103(zt,h0,h1,ione)
        call sk0(rhoj,dr,sk01)

        laph0 = -rhoj*rhoj*h0
        val = ejs(1)*rts(1)**3*(sk01-2*ima*h0)

        rhoj = -rts(2)
        call sk0(rhoj,dr,sk02)

        val = val - ejs(2)*rts(2)**3*(sk02)

        rhoj = -rts(3)
        call sk0(rhoj,dr,sk03)

        val = val - ejs(3)*rts(3)**3*(sk03)
        val = val/4

        alpha = 1/(ejs(1)*rts(1)**2+ejs(2)*rts(2)**2+ejs(3)*rts(3)**2)

        val = val - 1/(2*pi*alpha*dr)

        return
        end
c
c
c
c

c
c
c
c
c
        subroutine gradgshelm(rts,ejs,dx,dy,gradx,grady)
c
c  returns both components of grad G_s at displacement (dx,dy)
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(3),ejs(3),val,zt,ima,rhoj,h0,h1
        complex *16 sk01x,sk01y,sk02x,sk02y,sk03x,sk03y
        complex *16 gradx,grady,h0x,h0y
        real *8 dx,dy,dr,dr2,pi
        logical ilow
        data ima /(0,1)/

        dr = sqrt(dx**2+dy**2)
       
        zt = rts(1)*dr
        rhoj = rts(1)

        ione = 1
        call hank103(zt,h0,h1,ione)
        call gradsk0(rhoj,dx,dy,sk01x,sk01y)

        h0x = -rhoj*dx/dr*h1
        h0y = -rhoj*dy/dr*h1

        gradx = ejs(1)*rts(1)**2*(-sk01x+2*ima*h0x)
        grady = ejs(1)*rts(1)**2*(-sk01y+2*ima*h0y)

        rhoj = -rts(2)
        call gradsk0(rhoj,dx,dy,sk02x,sk02y)

        gradx = gradx + ejs(2)*rts(2)**2*sk02x
        grady = grady + ejs(2)*rts(2)**2*sk02y

        rhoj = -rts(3)
        call gradsk0(rhoj,dx,dy,sk03x,sk03y)

        gradx = gradx + ejs(3)*rts(3)**2*sk03x
        grady = grady + ejs(3)*rts(3)**2*sk03y
        
        gradx = gradx/2
        grady = grady/2

        return
        end
c
c
c
        subroutine gradxgshelm(rts,ejs,dx,dy,gradx)
c
c  returns the x-component of grad G_s at displacement (dx,dy)
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(3),ejs(3),val,zt,ima,rhoj,h0,h1
        complex *16 sk01x,sk01y,sk02x,sk02y,sk03x,sk03y
        complex *16 gradx,h0x,h0y
        real *8 dx,dy,dr,dr2,pi
        logical ilow
        data ima /(0,1)/

        dr = sqrt(dx**2+dy**2)

        zt = rts(1)*dr
        rhoj = rts(1)

        ione = 1
        call hank103(zt,h0,h1,ione)
        call gradsk0(rhoj,dx,dy,sk01x,sk01y)

        h0x = -rhoj*dx/dr*h1

        gradx = ejs(1)*rts(1)**2*(-sk01x+2*ima*h0x)

        rhoj = -rts(2)
        call gradsk0(rhoj,dx,dy,sk02x,sk02y)

        gradx = gradx + ejs(2)*rts(2)**2*sk02x

        rhoj = -rts(3)
        call gradsk0(rhoj,dx,dy,sk03x,sk03y)

        gradx = gradx + ejs(3)*rts(3)**2*sk03x

        gradx = gradx/2

        return
        end
c
c
c
c
c
        subroutine gradxgphihelm(rts,ejs,dx,dy,gradx)
c
c  returns the x-component of grad G_phi at displacement (dx,dy)
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(3),ejs(3),val,zt,ima,rhoj,h0,h1
        complex *16 sk01x,sk01y,sk02x,sk02y,sk03x,sk03y
        complex *16 gradx,h0x,h0y
        real *8 dx,dy,dr,dr2,pi
        logical ilow
        data ima /(0,1)/

        dr = sqrt(dx**2+dy**2)

        zt = rts(1)*dr
        rhoj = rts(1)

        ione = 1
        call hank103(zt,h0,h1,ione)
        call gradsk0(rhoj,dx,dy,sk01x,sk01y)

        h0x = -rhoj*dx/dr*h1

        gradx = ejs(1)*rts(1)**1*(-sk01x+2*ima*h0x)

        rhoj = -rts(2)
        call gradsk0(rhoj,dx,dy,sk02x,sk02y)

        gradx = gradx + ejs(2)*rts(2)**1*sk02x

        rhoj = -rts(3)
        call gradsk0(rhoj,dx,dy,sk03x,sk03y)

        gradx = gradx + ejs(3)*rts(3)**1*sk03x

        gradx = gradx/4

        return
        end
c
c
c
c
c
c
        subroutine gradygshelm(rts,ejs,dx,dy,grady)
c
c  returns the y-component of grad G_s at displacement (dx,dy)
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(3),ejs(3),val,zt,ima,rhoj,h0,h1
        complex *16 sk01x,sk01y,sk02x,sk02y,sk03x,sk03y
        complex *16 grady,h0x,h0y
        real *8 dx,dy,dr,dr2,pi
        logical ilow
        data ima /(0,1)/

        dr = sqrt(dx**2+dy**2)

        zt = rts(1)*dr
        rhoj = rts(1)

        ione = 1
        call hank103(zt,h0,h1,ione)
        call gradsk0(rhoj,dx,dy,sk01x,sk01y)

        h0y = -rhoj*dy/dr*h1

        grady = ejs(1)*rts(1)**2*(-sk01y+2*ima*h0y)

        rhoj = -rts(2)
        call gradsk0(rhoj,dx,dy,sk02x,sk02y)

        grady = grady + ejs(2)*rts(2)**2*sk02y

        rhoj = -rts(3)
        call gradsk0(rhoj,dx,dy,sk03x,sk03y)

        grady = grady + ejs(3)*rts(3)**2*sk03y

        grady = grady/2

        return
        end
c
c
c
c
c
c
        subroutine gradygphihelm(rts,ejs,dx,dy,grady)
c
c  returns the y-component of grad G_phi at displacement (dx,dy)
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(3),ejs(3),val,zt,ima,rhoj,h0,h1
        complex *16 sk01x,sk01y,sk02x,sk02y,sk03x,sk03y
        complex *16 grady,h0x,h0y
        real *8 dx,dy,dr,dr2,pi
        logical ilow
        data ima /(0,1)/

        dr = sqrt(dx**2+dy**2)

        zt = rts(1)*dr
        rhoj = rts(1)

        ione = 1
        call hank103(zt,h0,h1,ione)
        call gradsk0(rhoj,dx,dy,sk01x,sk01y)

        h0y = -rhoj*dy/dr*h1

        grady = ejs(1)*rts(1)**1*(-sk01y+2*ima*h0y)

        rhoj = -rts(2)
        call gradsk0(rhoj,dx,dy,sk02x,sk02y)

        grady = grady + ejs(2)*rts(2)**1*sk02y

        rhoj = -rts(3)
        call gradsk0(rhoj,dx,dy,sk03x,sk03y)

        grady = grady + ejs(3)*rts(3)**1*sk03y

        grady = grady/4

        return
        end
c
c
c        
c
c
c
        subroutine gradgphihelm(rts,ejs,dx,dy,gradx,grady)
c
c  returns both components of grad G_phi at displacement (dx,dy)
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(3),ejs(3),val,zt,ima,rhoj,h0,h1
        complex *16 sk01x,sk01y,sk02x,sk02y,sk03x,sk03y
        complex *16 gradx,grady,h0x,h0y
        real *8 dx,dy,dr,dr2,pi
        logical ilow
        data ima /(0,1)/

        dr = sqrt(dx**2+dy**2)

        zt = rts(1)*dr
        rhoj = rts(1)

        ione = 1
        call hank103(zt,h0,h1,ione)
        call gradsk0(rhoj,dx,dy,sk01x,sk01y)

        h0x = -rhoj*dx/dr*h1
        h0y = -rhoj*dy/dr*h1

        gradx = ejs(1)*rts(1)*(-sk01x+2*ima*h0x)
        grady = ejs(1)*rts(1)*(-sk01y+2*ima*h0y)

        rhoj = -rts(2)
        call gradsk0(rhoj,dx,dy,sk02x,sk02y)

        gradx = gradx + ejs(2)*rts(2)*sk02x
        grady = grady + ejs(2)*rts(2)*sk02y

        rhoj = -rts(3)
        call gradsk0(rhoj,dx,dy,sk03x,sk03y)

        gradx = gradx + ejs(3)*rts(3)*sk03x
        grady = grady + ejs(3)*rts(3)*sk03y

        gradx = gradx/4
        grady = grady/4

        return
        end
c
c
c
c
c
