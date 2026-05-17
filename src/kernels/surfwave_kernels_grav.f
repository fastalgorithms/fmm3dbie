!
! the following routines rely on the srcinfo and targinfo arrays
! containing the following information, standardized in the following
! way:
!
!     *info(1:3) = xyz
!     *info(4:5) = normal vector in the 2d plane
!
!     *zpars(1:3) = roots of cubic polynomial, 
!           zpars(1) must be the root closest to the real axis
!     *zpars(4:6) = residues
c
c
c
        subroutine gsgravkern(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,
     1    ipars,val)  
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
c
c
        subroutine gphigravkern(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,
     1    ipars,val)
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
        subroutine gsgrav(rts,ejs,dr,val)
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
c TODO ADD 1/r term
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
