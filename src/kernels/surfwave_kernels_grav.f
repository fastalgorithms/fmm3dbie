!
!  Gravity surface wave kernels. The *kern routines take zpars(1) = g,
!  with dispersion root rho = g/2 and residue 1. All values exclude the
!  algebraic term ej*rho/(2 pi r), which the near quadrature adds as a
!  scaled Laplace single layer.
!
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

        rts = zpars(1)/2.0d0
        ejs = 1.0d0

        call gsgrav(rts,ejs,dr,val)

        return
        end
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

        rts = zpars(1)/2.0d0
        ejs = 1.0d0

        call gphigrav(rts,ejs,dr,val)

        return
        end

c
        subroutine gsgrav(rts,ejs,dr,val)
c       G_S without the algebraic term:
c           (ej*rho^2/4) [ -K0(rho r) + 2i H0^(1)(rho r) ]
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(1),ejs(1),val,zt
        complex *16 h0,h1,sk01,ima,rhoj
        real *8 dr,pi
        data ima /(0,1)/

        rhoj = rts(1)
        val = 0

        if (abs(atan(imag(rhoj) / real(rhoj))) .le. 1d-14) then
c            real-positive (propagating) pole: sk0 and h0 = H0^(1)
             zt = rts(1)*dr
             ione = 1
             call hank103(zt,h0,h1,ione)
             call sk0(rhoj,dr,sk01)
             val = ejs(1)*rts(1)**2*(-sk01+2*ima*h0)
        else
c            complex/evanescent pole: sk0 only, no outgoing H0^(1) piece
             call sk0(-rhoj,dr,sk01)
             val = ejs(1)*rts(1)**2*sk01
        endif

        val = val/4.0d0


        return
        end

c
c
c
        subroutine gsgravgrad(rts,ejs,dx,dy,valx,valy)
c       Target gradient of gsgrav
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(1),ejs(1),valx,valy,ima,rhoj
        complex *16 sk0x,sk0y,h0,h1,zt,h0x,h0y
        real *8 dx,dy,dr
        data ima /(0,1)/

        rhoj = rts(1)
        dr = sqrt(dx**2+dy**2)

        if (abs(atan(imag(rhoj) / real(rhoj))) .le. 1d-14) then
c            real-positive (propagating) pole
             call gradsk0(rhoj,dx,dy,sk0x,sk0y)

             zt = rhoj*dr
             ione = 1
             call hank103(zt,h0,h1,ione)
c            d/dx H0^(1)(rho*r) = -rho*(dx/r)*H1^(1)(rho*r), likewise for y
             h0x = -rhoj*dx/dr*h1
             h0y = -rhoj*dy/dr*h1

             valx = ejs(1)*rts(1)**2*(-sk0x+2*ima*h0x)
             valy = ejs(1)*rts(1)**2*(-sk0y+2*ima*h0y)
        else
c            complex/evanescent pole: sk0(-rho,r) only
             call gradsk0(-rhoj,dx,dy,sk0x,sk0y)
             valx = ejs(1)*rts(1)**2*sk0x
             valy = ejs(1)*rts(1)**2*sk0y
        endif

        valx = valx/4.0d0
        valy = valy/4.0d0

        return
        end
c
c
c
        subroutine gphigravgrad(rts,ejs,dx,dy,valx,valy)
c       Target gradient of gphigrav
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(1),ejs(1),valx,valy,g

        call gsgravgrad(rts,ejs,dx,dy,valx,valy)
        g = 2.0d0*rts(1)
        valx = valx/g
        valy = valy/g

        return
        end
c
c
c
        subroutine gsgravgradxkern(src,ndt,targ,ndd,dpars,ndz,zpars,
     1    ndi,ipars,val)
c       x-component of the target gradient of gsgrav
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz
        complex *16 zpars(ndz), val, valx, valy
        complex *16 rts(1), ejs(1)

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)

        rts = zpars(1)/2.0d0
        ejs = 1.0d0

        call gsgravgrad(rts,ejs,dx,dy,valx,valy)
        val = valx

        return
        end
c
c
c
        subroutine gsgravgradykern(src,ndt,targ,ndd,dpars,ndz,zpars,
     1    ndi,ipars,val)
c       y-component of the target gradient of gsgrav
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz
        complex *16 zpars(ndz), val, valx, valy
        complex *16 rts(1), ejs(1)

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)

        rts = zpars(1)/2.0d0
        ejs = 1.0d0

        call gsgravgrad(rts,ejs,dx,dy,valx,valy)
        val = valy

        return
        end
c
c
c
        subroutine gphigravgradxkern(src,ndt,targ,ndd,dpars,ndz,
     1    zpars,ndi,ipars,val)
c       x-component of the target gradient of gphigrav
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz
        complex *16 zpars(ndz), val, valx, valy
        complex *16 rts(1), ejs(1)

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)

        rts = zpars(1)/2.0d0
        ejs = 1.0d0

        call gphigravgrad(rts,ejs,dx,dy,valx,valy)
        val = valx

        return
        end
c
c
c
        subroutine gphigravgradykern(src,ndt,targ,ndd,dpars,ndz,
     1    zpars,ndi,ipars,val)
c       y-component of the target gradient of gphigrav
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz
        complex *16 zpars(ndz), val, valx, valy
        complex *16 rts(1), ejs(1)

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)

        rts = zpars(1)/2.0d0
        ejs = 1.0d0

        call gphigravgrad(rts,ejs,dx,dy,valx,valy)
        val = valy

        return
        end
c
c


        subroutine gphigrav(rts,ejs,dr,val)
c
c  G_phi = G_S / g with g = 2*rho (see gsgrav)
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(1),ejs(1),val
        real *8 dr

        call gsgrav(rts,ejs,dr,val)
        val = val/(2*rts(1))

        return
        end
