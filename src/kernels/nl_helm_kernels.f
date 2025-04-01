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
!     *zpars(7:9) = alpha,beta,gamma in the exterior
!     *zpars(10:11) = beta,gamma in the interior
c
c
c
        subroutine gskern(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,
     1    ipars,val)  
        implicit real *8 (a-h,o-z)
        integer ndz, ndt
        real *8 src(*), targ(ndt)
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
c
c
c
c
        subroutine gshelm(rts,ejs,dr,val)
        implicit real *8 (a-h,o-z)
        complex *16 rts(3),ejs(3),val,zt
        complex *16 h0,h1,sk01,sk02,sk03,ima,rhoj
        real *8 dr,pi
        logical ilow
        data ima /(0,1)/

        zt = rts(1)*dr

        rhoj = rts(1)
        
        call hank103(zt,h0,h1,1)
        call sk0(rhoj,dr,sk01)
        
        val = ejs(1)*rts(1)**2*(-sk01+2*ima*h0)

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
c
        subroutine gphihelm(rts,ejs,dr,val)
        implicit real *8 (a-h,o-z)
        complex *16 rts(3),ejs(3),val,zt
        complex *16 h0,h1,sk01,sk02,sk03,ima,rhoj
        real *8 dr,pi
        logical ilow
        data ima /(0,1)/

        zt = rts(1)*dr

        rhoj = rts(1)

        call hank103(zt,h0,h1,1)
        call sk0(rhoj,dr,sk01)

        val = ejs(1)*rts(1)*(-sk01+2*ima*h0)

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
        subroutine lapgshelm(rts,ejs,dr,val)
        implicit real *8 (a-h,o-z)
        complex *16 rts(3),ejs(3),val,zt,laph0
        complex *16 h0,h1,lap1,lap2,lap3,ima,rhoj
        real *8 dr,pi
        logical ilow
        data ima /(0,1)/

        zt = rts(1)*dr

        rhoj = rts(1)

        call hank103(zt,h0,h1,1)
        call lapsk0(rhoj,dr,lap1)

        laph0 = -rhoj*rhoj*h0
        val = ejs(1)*rts(1)**2*(-lap1+2*ima*laph0)

        rhoj = -rts(2)
        call lapsk0(rhoj,dr,lap2)

        val = val + ejs(2)*rts(2)**2*lap2

        rhoj = -rts(3)
        call lapsk0(rhoj,dr,lap3)

        val = val + ejs(3)*rts(3)**2*lap3
        val = val/2

        return
        end
c
c
c
c
c
        subroutine lapgphihelm(rts,ejs,dr,val)
        implicit real *8 (a-h,o-z)
        complex *16 rts(3),ejs(3),val,zt,laph0
        complex *16 h0,h1,lap1,lap2,lap3,ima,rhoj
        real *8 dr,pi
        logical ilow  
        data ima /(0,1)/

        zt = rts(1)*dr

        rhoj = rts(1)

        call hank103(zt,h0,h1,1)
        call lapsk0(rhoj,dr,lap1)

        laph0 = -rhoj*rhoj*h0
        val = ejs(1)*rts(1)*(-lap1+2*ima*laph0)

        rhoj = -rts(2)
        call lapsk0(rhoj,dr,lap2)

        val = val + ejs(2)*rts(2)*lap2

        rhoj = -rts(3)   
        call lapsk0(rhoj,dr,lap3)

        val = val + ejs(3)*rts(3)*lap3
        val = val/4

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
        subroutine gradgshelm(rts,ejs,dx,dy,gradx,grady)
        implicit real *8 (a-h,o-z)
        complex *16 rts(3),ejs(3),val,zt,ima,rhoj,h0,h1
        complex *16 sk01x,sk01y,sk02x,sk02y,sk03x,sk03y
        complex *16 gradx,grady,h0x,h0y
        real *8 dx,dy,dr,dr2,pi
        logical ilow
        data ima /(0,1)/

        dr = sqrt(dx**2+dy**2)
       
        zt = rts(1)*dr
        rhoj = rts(1)

        call hank103(zt,h0,h1,1)
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
c
c
c
        subroutine gradgphihelm(rts,ejs,dx,dy,gradx,grady)
        implicit real *8 (a-h,o-z)
        complex *16 rts(3),ejs(3),val,zt,ima,rhoj,h0,h1
        complex *16 sk01x,sk01y,sk02x,sk02y,sk03x,sk03y
        complex *16 gradx,grady,h0x,h0y
        real *8 dx,dy,dr,dr2,pi
        logical ilow
        data ima /(0,1)/

        dr = sqrt(dx**2+dy**2)

        zt = rts(1)*dr
        rhoj = rts(1)

        call hank103(zt,h0,h1,1)
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
