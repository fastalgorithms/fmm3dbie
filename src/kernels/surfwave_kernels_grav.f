! Kernels for gravity surface wave problems
!
! the following routines rely on the srcinfo and targinfo arrays
! containing the following information, standardized in the following
! way:
!
!     *info(1:3) = xyz
!     *info(4:5) = tangent vector in the 2d plane
!     *info(10:11) = normal vector in the 2d plane
!
!     *zpars(1) = root of the dispersion relation
!     *zpars(2) = residue (set to 1.0 in kernel wrappers)
!
! The gravity Green's functions gs and gphi are evaluated via
! the primitive routines gsgrav and gphigrav.
!
c
c
c
c  Kernel wrappers (src/targ in srcinfo/targinfo format)
c
c
        subroutine gsgravkern(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,
     1    ipars,val)
c
c  returns the gravity wave free-surface Green's function G_s(r)
c  evaluated at r = |x-y|, using zpars(1) as the dispersion root
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz,dr
        complex *16 zpars(ndz), val        
        complex *16 rts(1), ejs(1)

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)

        dr=sqrt(dx**2+dy**2+dz**2)

        rts = zpars(1)
        ejs = 1.0d0

        call gsgrav(rts,ejs,dr,val)

        return 
        end
c
c
c
        subroutine gphigravkern(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,
     1    ipars,val)
c
c  returns the gravity wave velocity potential Green's function G_phi(r)
c  evaluated at r = |x-y|, using zpars(1) as the dispersion root
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz,dr
        complex *16 zpars(ndz), val
        complex *16 rts(1), ejs(1)

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)

        dr=sqrt(dx**2+dy**2+dz**2)

        rts = zpars(1)
        ejs = 1.0d0

        call gphigrav(rts,ejs,dr,val)

        return
        end

c
c
c
c
c        subroutine s3dgphigravkern(src,ndt,targ,ndd,dpars,ndz,zpars,
c     1    ndi,ipars,val)
c        implicit real *8 (a-h,o-z)
c        integer ndz, ndt
c        real *8 src(*), targ(ndt)
c        real *8 dx,dy,dz,dr
c        complex *16 zpars(ndz), val
c        complex *16 rts(3), ejs(3)
c
c        dx=targ(1)-src(1)
c        dy=targ(2)-src(2)
c        dz=targ(3)-src(3)
c
c        dr=sqrt(dx**2+dy**2+dz**2)
c
c        rts = 1.0d0
c        ejs = 1.0d0 
c
c        call s3dgphigrav(rts,ejs,dr,val)
c
c        return
c        end
c
c
c
c
c
c
c
c
c
c  Primitive Green's function evaluators (r-only interface)
c
        subroutine gsgrav(rts,ejs,dr,val)
c
c  returns the gravity wave free-surface Green's function G_s at distance dr
c    G_s(r) = ejs(1)*rts(1)^2 * (-K_0(rts(1)*r) + 2i*H_0(rts(1)*r)) / 2
c  where rts(1) is the dispersion root on the real axis
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(1),ejs(1),val,zt
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
        
        val = val/2.0d0
        
c        val = val + 1.0d0/dr/pi

        return
        end
c
c
c
c
        subroutine gphigrav(rts,ejs,dr,val)
c
c  returns the gravity wave velocity potential Green's function G_phi at
c  distance dr. For the real dispersion root uses Hankel/K_0; for the
c  imaginary root uses K_0 with negated argument.
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(1),ejs(1),val,zt
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

        val = val/4
c        val = val + 0.5d0/dr/pi

        return
        end
c
c
c
        subroutine s3dgphigrav(rts,ejs,dr,val)
c
c  returns the 3D free-space contribution to the gravity wave velocity
c  potential Green's function at distance dr (same form as gphigrav
c  but without the rts(1) prefactor, scaled by 1/8)
c  TODO: add 1/r free-space term
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(1),ejs(1),val,zt
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

        val = val/8

        return
        end
c
c
