!
!
!
! the following routine is used to generate derivatives of the 
! Hankel function for various source target pairs
!
!
        subroutine helmdiffgreen(zk,dx,dy,h0,h0x,h0y,h0xx,
     1      h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy)
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        real *8 dx,dy,dx2,dy2,dx3,dy3,dx4,dy4
        real *8 dr,r2,r,rm1,rm2,rm3,rm4
        integer *8 nder
        complex *16 d0,d1,d2,d3,ima
        complex *16 h0,h0x,h0y,zk,zt,zk2,zk3
        complex *16 h0xx,h0xy,h0yy
        complex *16 h0xxx,h0xxy,h0xyy,h0yyy
        complex *16 cvals(4,1)
        data ima /(0,1)/
        real *8 pi

        dr = sqrt(dx**2+dy**2)

        dx2 = dx*dx
        dx3 = dx2*dx
        dx4 = dx3*dx

        dy2 = dy*dy
        dy3 = dy2*dy
        dy4 = dy3*dy

        rm1 = 1/dr
        rm2 = rm1*rm1
        rm3 = rm2*rm1
        rm4 = rm3*rm1
        rm5 = rm4*rm1

        zk2 = zk*zk
        zk3 = zk2*zk


        zt = zk*dr
        ione = 1
        ithree = 3
        call hankdiff(zt,ione,ithree,ione,cvals)

        d0 = cvals(1,1)
        d1 = cvals(2,1)
        d2 = cvals(3,1) 
        d3 = cvals(4,1)

        h0x = zk*dx*d1*rm1
        h0y = zk*dy*d1*rm1

        h0xx = zk*dy2*rm3*d1 + zk2*dx2*rm2*d2 
        h0xy = -zk*dx*dy*rm3*d1 + zk2*dx*dy*rm2*d2 
        h0yy = zk*dx2*rm3*d1 + zk2*dy2*rm2*d2

        h0xxx = (-3*zk*dy2*dx*rm5*d1 + zk2*dx*dy2*rm4*d2  
     1             + 2*zk2*dy2*dx*rm4*d2 + zk3*dx3*rm3*d3) 
        h0xxy = (-zk*(dy3 - 2*dx2*dy)*rm5*d1 + zk2*dy3*rm4*d2  
     1             - 2*zk2*dx2*dy*rm4*d2 + zk3*dx2*dy*rm3*d3) 
        h0xyy = (-zk*(dx3 - 2*dy2*dx)*rm5*d1 + zk2*dx3*rm4*d2 
     1             - 2*zk2*dy2*dx*rm4*d2 + zk3*dy2*dx*rm3*d3)
        h0yyy = (-3*zk*dx2*dy*rm5*d1 + zk2*dy*dx2*rm4*d2  
     1             + 2*zk2*dx2*dy*rm4*d2 + zk3*dy3*rm3*d3) 
        
        pi = atan(1d0)*4
        h0 = ima*d0/4-log(zk)/(2*pi)
c        h0 = ima*h0/4

        h0x = ima*h0x/4
        h0y = ima*h0y/4

        h0xx = ima*h0xx/4
        h0xy = ima*h0xy/4
        h0yy = ima*h0yy/4

        h0xxx = ima*h0xxx/4
        h0xxy = ima*h0xxy/4
        h0xyy = ima*h0xyy/4
        h0yyy = ima*h0yyy/4

        return
        end
c
       subroutine helmdiffgreen23(zk,dx,dy,h0xx,
     1      h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy)
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        real *8 dx,dy,dx2,dy2,dx3,dy3,dx4,dy4
        real *8 dr,r2,r,rm1,rm2,rm3,rm4
        integer *8 nder
        complex *16 d0,d1,d2,d3,ima
        complex *16 h0,h0x,h0y,zk,zt,zk2,zk3
        complex *16 h0xx,h0xy,h0yy
        complex *16 h0xxx,h0xxy,h0xyy,h0yyy
        complex *16 cvals(4,1)
        data ima /(0,1)/

        dr = sqrt(dx**2+dy**2)

        dx2 = dx*dx
        dx3 = dx2*dx
        dx4 = dx3*dx

        dy2 = dy*dy
        dy3 = dy2*dy
        dy4 = dy3*dy

        rm1 = 1/dr
        rm2 = rm1*rm1
        rm3 = rm2*rm1
        rm4 = rm3*rm1
        rm5 = rm4*rm1

        zk2 = zk*zk
        zk3 = zk2*zk


        zt = zk*dr
        ione = 1
        ithree = 3
        call hankdiff(zt,ione,ithree,ione,cvals)




        d0 = cvals(1,1)
        d1 = cvals(2,1)
        d2 = cvals(3,1) 
        d3 = cvals(4,1)


        h0xx = zk*dy2*rm3*d1 + zk2*dx2*rm2*d2 
        h0xy = -zk*dx*dy*rm3*d1 + zk2*dx*dy*rm2*d2 
        h0yy = zk*dx2*rm3*d1 + zk2*dy2*rm2*d2

        h0xxx = (-3*zk*dy2*dx*rm5*d1 + zk2*dx*dy2*rm4*d2  
     1             + 2*zk2*dy2*dx*rm4*d2 + zk3*dx3*rm3*d3) 
        h0xxy = (-zk*(dy3 - 2*dx2*dy)*rm5*d1 + zk2*dy3*rm4*d2  
     1             - 2*zk2*dx2*dy*rm4*d2 + zk3*dx2*dy*rm3*d3) 
        h0xyy = (-zk*(dx3 - 2*dy2*dx)*rm5*d1 + zk2*dx3*rm4*d2 
     1             - 2*zk2*dy2*dx*rm4*d2 + zk3*dy2*dx*rm3*d3)
        h0yyy = (-3*zk*dx2*dy*rm5*d1 + zk2*dy*dx2*rm4*d2  
     1             + 2*zk2*dx2*dy*rm4*d2 + zk3*dy3*rm3*d3) 
        
        h0xx = ima*h0xx/4
        h0xy = ima*h0xy/4
        h0yy = ima*h0yy/4

        h0xxx = ima*h0xxx/4
        h0xxy = ima*h0xxy/4
        h0xyy = ima*h0xyy/4
        h0yyy = ima*h0yyy/4

        return
        end
c
c
c
c
        subroutine getders(zt,d0,d1,d2,d3)
        complex *16 zt,h0,h1
        complex *16 d0,d1,d2,d3

        ione = 1
        call hank103(zt,h0,h1,ione)

        d0 = h0
        d1 = -h1
        d2 = -h0+h1/zt
        d3 = h1- h1/(zt**2) + (h0 - 1/zt*h1)/zt
        
        return
        end
c
c
