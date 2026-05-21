! Kernels for flexural surface wave problems
!
! the following routines rely on the srcinfo and targinfo arrays
! containing the following information, standardized in the following
! way:
!
!     *info(1:3) = xyz
!     *info(4:5) = tangent vector in the 2d plane
!     *info(10:11) = normal vector in the 2d plane
!     *info(13) = kappa (curvature)
!
!     *zpars(1:5) = roots of the quintic dispersion polynomial,
!           zpars(1) must be the root closest to the real axis
!     *zpars(6:10) = residues corresponding to zpars(1:5)
!     *zpars(11) = alpha (stiffness parameter)
!     *zpars(12) = gamma (mass loading parameter)
!     *zpars(13) = nu (Poisson ratio)
!
! The flexural Green's functions gs and gphi are weighted sums over
! the five dispersion roots, evaluated via gsflex / gphiflex (r-only)
! or gsflex23 / gphiflex23 (second and third derivative tensors).
!
c
c
c  Kernel wrappers (src/targ in srcinfo/targinfo format)
c
c
        subroutine flexrepbc1(srcinfo,ndt,targinfo,ndd,dpars,ndz,
     1    zpars,ndi,ipars,val)
c
c  returns the first boundary condition kernel for the flexural rep:
c  moment condition (2nd-order operator on gs and gphi weighted by gamma and nu)
c    val = (-gamma/2 * M[gs] + M[gphi]) / 2
c  where M[g] = (n^T Hess(g) n) + nu*(tau^T Hess(g) tau)
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 srcinfo(*), targinfo(ndt)
        real *8 dx,dy,dz,dr,nx,ny,taux,tauy
        complex *16 zpars(ndz), val, val1, val2
        complex *16 rts(5), ejs(5),alpha,beta,gamma,nu
        complex *16 gsxx,gsxy,gsyy,gsxxx,gsxxy,gsxyy,gsyyy,gphiyyy
        complex *16 gphixx,gphixy,gphiyy,gphixxx,gphixxy,gphixyy

        dx=targinfo(1)-srcinfo(1)
        dy=targinfo(2)-srcinfo(2)

        dr=sqrt(dx**2+dy**2)

        taux = targinfo(4)
        tauy = targinfo(5)

        nx = targinfo(10)
        ny = targinfo(11)

        rts = zpars(1:5)
        ejs = zpars(6:10)

        alpha = zpars(11)
        gamma = zpars(12)
        nu = zpars(13)

        call gsflex23(rts,ejs,dx,dy,gsxx,gsxy,gsyy,gsxxx,gsxxy,
     1     gsxyy,gsyyy)

        call gphiflex23(rts,ejs,dx,dy,gphixx,gphixy,gphiyy,gphixxx,
     1     gphixxy,gphixyy,gphiyyy)

        val1 = gsxx*nx*nx + 2*gsxy*nx*ny + gsyy*ny*ny
     1       + nu*(gsxx*taux*taux+2*gsxy*taux*tauy+gsyy*tauy*tauy)

        val2 = gphixx*nx*nx + 2*gphixy*nx*ny + gphiyy*ny*ny
     1     + nu*(gphixx*taux*taux+2*gphixy*taux*tauy+gphiyy*tauy*tauy)

        val = -gamma/2*val1 + val2
        val = val/2

        return
        end
c
c
c
c
        subroutine flexrepbc2(srcinfo,ndt,targinfo,ndd,dpars,ndz,
     1    zpars,ndi,ipars,val)
c
c  returns the second boundary condition kernel for the flexural rep:
c  shear force condition (3rd-order operator on gs and gphi with kappa)
c    val = (-gamma/2 * V[gs] + V[gphi]) / 2
c  where V[g] is the Kirchhoff shear force involving kappa, n, tau, and nu
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 srcinfo(*), targinfo(ndt)
        real *8 dx,dy,dz,dr,nx,ny,taux,tauy,kappa
        complex *16 zpars(ndz), val, val1, val2
        complex *16 rts(5), ejs(5),alpha,beta,gamma,nu
        complex *16 gsxx,gsxy,gsyy,gsxxx,gsxxy,gsxyy,gsyyy,gphiyyy
        complex *16 gphixx,gphixy,gphiyy,gphixxx,gphixxy,gphixyy

        dx=targinfo(1)-srcinfo(1)
        dy=targinfo(2)-srcinfo(2)

        dr=sqrt(dx**2+dy**2)

        taux = targinfo(4)
        tauy = targinfo(5)

        nx = targinfo(10)
        ny = targinfo(11)

        kappa = targinfo(13)
        
        rts = zpars(1:5)
        ejs = zpars(6:10)

        alpha = zpars(11)
        gamma = zpars(12)
        nu = zpars(13)

        call gsflex23(rts,ejs,dx,dy,gsxx,gsxy,gsyy,gsxxx,gsxxy,
     1     gsxyy,gsyyy)

        call gphiflex23(rts,ejs,dx,dy,gphixx,gphixy,gphiyy,gphixxx,
     1     gphixxy,gphixyy,gphiyyy)

        val1 = gsxxx*nx*nx*nx + 3*gsxxy*nx*nx*ny + 3*gsxyy*nx*ny*ny
     1      + gsyyy*ny*ny*ny + (2-nu)*(gsxxx*nx*taux*taux
     1      + gsxxy*(taux*taux*ny + 2*taux*tauy*nx)
     1      + gsxyy*(2*taux*tauy*ny + tauy*tauy*nx)
     1      + gsyyy*tauy*tauy*ny)+kappa*(1-nu)*(gsxx*taux*taux
     1      + 2*gsxy*taux*tauy+gsyy*tauy*tauy - gsxx*nx*nx
     1      - 2*gsxy*nx*ny - gsyy*ny*ny)

        val2 = gphixxx*nx*nx*nx + 3*gphixxy*nx*nx*ny
     1      + 3*gphixyy*nx*ny*ny
     1      + gphiyyy*ny*ny*ny + (2-nu)*(gphixxx*nx*taux*taux
     1      + gphixxy*(taux*taux*ny + 2*taux*tauy*nx)
     1      + gphixyy*(2*taux*tauy*ny + tauy*tauy*nx)
     1      + gphiyyy*tauy*tauy*ny)+kappa*(1-nu)*(gphixx*taux*taux
     1      + 2*gphixy*taux*tauy+gphiyy*tauy*tauy - gphixx*nx*nx
     1      - 2*gphixy*nx*ny - gphiyy*ny*ny)

        val = -gamma/2*val1 + val2
        val = val/2

        return
        end
c
c
c
        subroutine gsflexkern(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,
     1    ipars,val)
c
c  returns the flexural free-surface Green's function G_s(r) at r = |x-y|
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz,dr
        complex *16 zpars(ndz), val        
        complex *16 rts(5), ejs(5)

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)

        dr=sqrt(dx**2+dy**2+dz**2)

        rts = zpars(1:5)
        ejs = zpars(6:10) 

        call gsflex(rts,ejs,dr,val)

        return 
        end
c
c
c
        subroutine bilapgsflexkern(src,ndt,targ,ndd,dpars,ndz,zpars,
     1    ndi,ipars,val)
c
c  returns the biharmonic (Delta^2) of the free-surface Green's function:
c  Delta^2 G_s(r)
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz,dr
        complex *16 zpars(ndz), val
        complex *16 rts(5), ejs(5)

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)

        dr=sqrt(dx**2+dy**2+dz**2)

        rts = zpars(1:5)
        ejs = zpars(6:10)

        call bilapgsflex(rts,ejs,dr,val)

        return
        end
c
c
c
        subroutine s3dgphiflexkern(src,ndt,targ,ndd,dpars,ndz,zpars,
     1    ndi,ipars,val)
c
c  returns the 3D free-space part of the velocity potential Green's
c  function G_phi(r): residue-weighted sum without rts prefactor, scaled by 1/8
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz,dr
        complex *16 zpars(ndz), val
        complex *16 rts(5), ejs(5)

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)

        dr=sqrt(dx**2+dy**2+dz**2)

        rts = zpars(1:5)
        ejs = zpars(6:10)

        call s3dgphiflex(rts,ejs,dr,val)

        return
        end
c
c
c        
c
c
        subroutine bilapgphiflexkern(src,ndt,targ,ndd,dpars,ndz,zpars,
     1    ndi,ipars,val)
c
c  returns the biharmonic of the velocity potential Green's function:
c  Delta^2 G_phi(r), including the 1/r free-space correction
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz,dr
        complex *16 zpars(ndz), val
        complex *16 rts(5), ejs(5)

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)

        dr=sqrt(dx**2+dy**2+dz**2)

        rts = zpars(1:5)
        ejs = zpars(6:10)

        call bilapgphiflex(rts,ejs,dr,val)

        return
        end
c
c
c
c
c
c
        subroutine gphiflexkern(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,
     1    ipars,val)
c
c  returns the flexural velocity potential Green's function G_phi(r)
c  at r = |x-y|
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        integer *8 ndz, ndt
        real *8 src(*), targ(ndt)
        real *8 dx,dy,dz,dr
        complex *16 zpars(ndz), val
        complex *16 rts(5), ejs(5)

        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)

        dr=sqrt(dx**2+dy**2+dz**2)

        rts = zpars(1:5)
        ejs = zpars(6:10)

        call gphiflex(rts,ejs,dr,val)

        return
        end

c
c
c
c
c  Primitive Green's function evaluators
c
        subroutine gsflex(rts,ejs,dr,val)
c
c  returns the flexural free-surface Green's function G_s at distance dr:
c  weighted sum over 5 roots, real root via Hankel+K_0, complex roots via K_0
c  G_s(r) = sum_j ejs(j)*rts(j)^2 * (-K_0 + 2i*H_0 or K_0) / 2
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(5),ejs(5),val,zt
        complex *16 h0,h1,sk01,sk02,sk03,sk04,sk05
        complex *16 ima,rhoj
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

        rhoj = -rts(4)
        call sk0(rhoj,dr,sk04)

        val = val + ejs(4)*rts(4)**2*sk04

        rhoj = -rts(5)
        call sk0(rhoj,dr,sk05)

        val = val + ejs(5)*rts(5)**2*sk05

        val = val/2

        return
        end
c
c
c
        subroutine gsflex23(rts,ejs,dx,dy,gsxx,gsxy,gsyy,gsxxx,gsxxy,
     1     gsxyy,gsyyy)
c
c  returns the 2nd and 3rd derivative tensors of G_s at displacement (dx,dy):
c  gsxx, gsxy, gsyy (Hessian) and gsxxx, gsxxy, gsxyy, gsyyy (3rd order)
c  uses helmdiffgreen23 for real root and struveKdiffgreen23 for all roots
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(5),ejs(5),zt,gsxx,gsxy,gsyy
        complex *16 gsxxx,gsxxy,gsxyy,gsyyy
        complex *16 h0xx,h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy
        complex *16 sk0xx,sk0xy,sk0yy,sk0xxx,sk0xxy,sk0xyy,sk0yyy
        complex *16 ima,rhoj
        real *8 dx,dy,dr,pi,dx2,dx3,dy2,dy3,dr2,dr3
        logical ilow
        data ima /(0,1)/

        zt = rts(1)*dr

        rhoj = rts(1)

        call helmdiffgreen23(rhoj,dx,dy,h0xx,
     1      h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy)

        h0xx = -4*ima*h0xx
        h0xy = -4*ima*h0xy
        h0yy = -4*ima*h0yy
        h0xxx = -4*ima*h0xxx
        h0xxy = -4*ima*h0xxy
        h0xyy = -4*ima*h0xyy
        h0yyy = -4*ima*h0yyy
    
        call struveKdiffgreen23(rhoj,dx,dy,sk0xx,sk0xy,sk0yy,
     1      sk0xxx,sk0xxy,sk0xyy,sk0yyy)

        gsxx = ejs(1)*rts(1)**2*(-sk0xx + 2*ima*h0xx)
        gsxy = ejs(1)*rts(1)**2*(-sk0xy + 2*ima*h0xy)
        gsyy = ejs(1)*rts(1)**2*(-sk0yy + 2*ima*h0yy)
    
        gsxxx = ejs(1)*rts(1)**2*(-sk0xxx + 2*ima*h0xxx)
        gsxxy = ejs(1)*rts(1)**2*(-sk0xxy + 2*ima*h0xxy)
        gsxyy = ejs(1)*rts(1)**2*(-sk0xyy + 2*ima*h0xyy)
        gsyyy = ejs(1)*rts(1)**2*(-sk0yyy + 2*ima*h0yyy)


        rhoj = -rts(2)
        call struveKdiffgreen23(rhoj,dx,dy,sk0xx,sk0xy,sk0yy,
     1      sk0xxx,sk0xxy,sk0xyy,sk0yyy)


        gsxx = gsxx + ejs(2)*rts(2)**2*(sk0xx)
        gsxy = gsxy + ejs(2)*rts(2)**2*(sk0xy)
        gsyy = gsyy + ejs(2)*rts(2)**2*(sk0yy)

        gsxxx = gsxxx + ejs(2)*rts(2)**2*(sk0xxx)
        gsxxy = gsxxy + ejs(2)*rts(2)**2*(sk0xxy)
        gsxyy = gsxyy + ejs(2)*rts(2)**2*(sk0xyy)
        gsyyy = gsyyy + ejs(2)*rts(2)**2*(sk0yyy)


        rhoj = -rts(3)

        call struveKdiffgreen23(rhoj,dx,dy,sk0xx,sk0xy,sk0yy,
     1      sk0xxx,sk0xxy,sk0xyy,sk0yyy)

        gsxx = gsxx + ejs(3)*rts(3)**2*(sk0xx)
        gsxy = gsxy + ejs(3)*rts(3)**2*(sk0xy)
        gsyy = gsyy + ejs(3)*rts(3)**2*(sk0yy)

        gsxxx = gsxxx + ejs(3)*rts(3)**2*(sk0xxx)
        gsxxy = gsxxy + ejs(3)*rts(3)**2*(sk0xxy)
        gsxyy = gsxyy + ejs(3)*rts(3)**2*(sk0xyy)
        gsyyy = gsyyy + ejs(3)*rts(3)**2*(sk0yyy)

        rhoj = -rts(4)

        call struveKdiffgreen23(rhoj,dx,dy,sk0xx,sk0xy,sk0yy,
     1      sk0xxx,sk0xxy,sk0xyy,sk0yyy)


        gsxx = gsxx + ejs(4)*rts(4)**2*(sk0xx)
        gsxy = gsxy + ejs(4)*rts(4)**2*(sk0xy)
        gsyy = gsyy + ejs(4)*rts(4)**2*(sk0yy)

        gsxxx = gsxxx + ejs(4)*rts(4)**2*(sk0xxx)
        gsxxy = gsxxy + ejs(4)*rts(4)**2*(sk0xxy)
        gsxyy = gsxyy + ejs(4)*rts(4)**2*(sk0xyy)
        gsyyy = gsyyy + ejs(4)*rts(4)**2*(sk0yyy)

        rhoj = -rts(5)

        call struveKdiffgreen23(rhoj,dx,dy,sk0xx,sk0xy,sk0yy,
     1      sk0xxx,sk0xxy,sk0xyy,sk0yyy)


        gsxx = gsxx + ejs(5)*rts(5)**2*(sk0xx)
        gsxy = gsxy + ejs(5)*rts(5)**2*(sk0xy)
        gsyy = gsyy + ejs(5)*rts(5)**2*(sk0yy)

        gsxxx = gsxxx + ejs(5)*rts(5)**2*(sk0xxx)
        gsxxy = gsxxy + ejs(5)*rts(5)**2*(sk0xxy)
        gsxyy = gsxyy + ejs(5)*rts(5)**2*(sk0xyy)
        gsyyy = gsyyy + ejs(5)*rts(5)**2*(sk0yyy)


        gsxx = gsxx/2
        gsxy = gsxy/2
        gsyy = gsyy/2
        gsxxx = gsxxx/2
        gsxxy = gsxxy/2
        gsxyy = gsxyy/2
        gsyyy = gsyyy/2

        return
        end
c
c
        subroutine gphiflex23(rts,ejs,dx,dy,gsxx,gsxy,gsyy,gsxxx,gsxxy,
     1     gsxyy,gsyyy)
c
c  returns the 2nd and 3rd derivative tensors of G_phi at displacement (dx,dy)
c  same structure as gsflex23 but with rts(j)^1 prefactor; includes a
c  1/(2*pi) correction term in the 3rd derivatives from the 1/r free-space term
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(5),ejs(5),zt,gsxx,gsxy,gsyy
        complex *16 gsxxx,gsxxy,gsxyy,gsyyy
        complex *16 h0xx,h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy
        complex *16 sk0xx,sk0xy,sk0yy,sk0xxx,sk0xxy,sk0xyy,sk0yyy
        complex *16 ima,rhoj,oneoveralpha
        real *8 dx,dy,dr,pi,dx2,dy2,dx3,dy3,dr2,dr3
        logical ilow
        data ima /(0,1)/,pi/0.31415926535897932D+01/

        oneoveralpha = ejs(1)*rts(1)**4+ejs(2)*rts(2)**4
     1           +ejs(3)*rts(3)**4+ejs(4)*rts(4)**4
     1           +ejs(5)*rts(5)**4

        dr = sqrt(dx**2+dy**2)

        zt = rts(1)*dr

        rhoj = rts(1)

        call helmdiffgreen23(rhoj,dx,dy,h0xx,
     1      h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy)

        h0xx = -4*ima*h0xx
        h0xy = -4*ima*h0xy
        h0yy = -4*ima*h0yy
        h0xxx = -4*ima*h0xxx
        h0xxy = -4*ima*h0xxy
        h0xyy = -4*ima*h0xyy
        h0yyy = -4*ima*h0yyy
    
        call struveKdiffgreen23(rhoj,dx,dy,sk0xx,sk0xy,sk0yy,
     1      sk0xxx,sk0xxy,sk0xyy,sk0yyy)

        gsxx = ejs(1)*rts(1)*(-sk0xx + 2*ima*h0xx)
        gsxy = ejs(1)*rts(1)*(-sk0xy + 2*ima*h0xy)
        gsyy = ejs(1)*rts(1)*(-sk0yy + 2*ima*h0yy)
    
        gsxxx = ejs(1)*rts(1)*(-sk0xxx + 2*ima*h0xxx)
        gsxxy = ejs(1)*rts(1)*(-sk0xxy + 2*ima*h0xxy)
        gsxyy = ejs(1)*rts(1)*(-sk0xyy + 2*ima*h0xyy)
        gsyyy = ejs(1)*rts(1)*(-sk0yyy + 2*ima*h0yyy)


        rhoj = -rts(2)
        call struveKdiffgreen23(rhoj,dx,dy,sk0xx,sk0xy,sk0yy,
     1      sk0xxx,sk0xxy,sk0xyy,sk0yyy)


        gsxx = gsxx + ejs(2)*rts(2)*(sk0xx)
        gsxy = gsxy + ejs(2)*rts(2)*(sk0xy)
        gsyy = gsyy + ejs(2)*rts(2)*(sk0yy)

        gsxxx = gsxxx + ejs(2)*rts(2)*(sk0xxx)
        gsxxy = gsxxy + ejs(2)*rts(2)*(sk0xxy)
        gsxyy = gsxyy + ejs(2)*rts(2)*(sk0xyy)
        gsyyy = gsyyy + ejs(2)*rts(2)*(sk0yyy)


        rhoj = -rts(3)

        call struveKdiffgreen23(rhoj,dx,dy,sk0xx,sk0xy,sk0yy,
     1      sk0xxx,sk0xxy,sk0xyy,sk0yyy)

        gsxx = gsxx + ejs(3)*rts(3)*(sk0xx)
        gsxy = gsxy + ejs(3)*rts(3)*(sk0xy)
        gsyy = gsyy + ejs(3)*rts(3)*(sk0yy)

        gsxxx = gsxxx + ejs(3)*rts(3)*(sk0xxx)
        gsxxy = gsxxy + ejs(3)*rts(3)*(sk0xxy)
        gsxyy = gsxyy + ejs(3)*rts(3)*(sk0xyy)
        gsyyy = gsyyy + ejs(3)*rts(3)*(sk0yyy)

        rhoj = -rts(4)

        call struveKdiffgreen23(rhoj,dx,dy,sk0xx,sk0xy,sk0yy,
     1      sk0xxx,sk0xxy,sk0xyy,sk0yyy)


        gsxx = gsxx + ejs(4)*rts(4)*(sk0xx)
        gsxy = gsxy + ejs(4)*rts(4)*(sk0xy)
        gsyy = gsyy + ejs(4)*rts(4)*(sk0yy)

        gsxxx = gsxxx + ejs(4)*rts(4)*(sk0xxx)
        gsxxy = gsxxy + ejs(4)*rts(4)*(sk0xxy)
        gsxyy = gsxyy + ejs(4)*rts(4)*(sk0xyy)
        gsyyy = gsyyy + ejs(4)*rts(4)*(sk0yyy)

        rhoj = -rts(5)

        call struveKdiffgreen23(rhoj,dx,dy,sk0xx,sk0xy,sk0yy,
     1      sk0xxx,sk0xxy,sk0xyy,sk0yyy)


        gsxx = gsxx + ejs(5)*rts(5)*(sk0xx)
        gsxy = gsxy + ejs(5)*rts(5)*(sk0xy)
        gsyy = gsyy + ejs(5)*rts(5)*(sk0yy)

        gsxxx = gsxxx + ejs(5)*rts(5)*(sk0xxx)
        gsxxy = gsxxy + ejs(5)*rts(5)*(sk0xxy)
        gsxyy = gsxyy + ejs(5)*rts(5)*(sk0xyy)
        gsyyy = gsyyy + ejs(5)*rts(5)*(sk0yyy)


        dx2 = dx*dx
        dx3 = dx2*dx

        dy2 = dy*dy
        dy3 = dy2*dy
        
        dr2 = dr*dr
        dr3 = dr2*dr

        gsxx = gsxx/4
        gsxy = gsxy/4
        gsyy = gsyy/4
        gsxxx = oneoveralpha/(2*pi)*dx3/dr3 + gsxxx/4
        gsxxy = oneoveralpha/(2*pi)*dx2*dy/dr3 + gsxxy/4
        gsxyy = oneoveralpha/(2*pi)*dx*dy2/dr3 + gsxyy/4
        gsyyy = oneoveralpha/(2*pi)*dy3/dr3 + gsyyy/4

        return
        end
c
c
c
c
        subroutine gphiflex(rts,ejs,dr,val)
c
c  returns the flexural velocity potential Green's function G_phi at distance dr:
c  same structure as gsflex but with rts(j)^1 prefactor, scaled by 1/4
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(5),ejs(5),val,zt
        complex *16 h0,h1,sk01,sk02,sk03,sk04,sk05
        complex *16 ima,rhoj
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

        rhoj = -rts(4)
        call sk0(rhoj,dr,sk04)

        val = val + ejs(4)*rts(4)*sk04

        rhoj = -rts(5)
        call sk0(rhoj,dr,sk05)

        val = val + ejs(5)*rts(5)*sk05

        val = val/4

        return
        end
c
c
c
c
c
        subroutine bilapgsflex(rts,ejs,dr,val)
c
c  returns Delta^2 G_s(r): biharmonic of the free-surface Green's function.
c  rts(j)^6 prefactor, real root uses Hankel+K_0, complex roots use K_0.
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(5),ejs(5),val,zt
        complex *16 h0,h1,sk01,sk02,sk03,sk04,sk05
        complex *16 ima,rhoj
        real *8 dr,pi
        logical ilow
        data ima /(0,1)/

        zt = rts(1)*dr

        rhoj = rts(1)
        
        ione = 1
        call hank103(zt,h0,h1,ione)
        call sk0(rhoj,dr,sk01)
        
        val = ejs(1)*rts(1)**6*(-sk01+2*ima*h0)

        rhoj = -rts(2)
        call sk0(rhoj,dr,sk02)

        val = val + ejs(2)*rts(2)**6*sk02
        
        rhoj = -rts(3)
        call sk0(rhoj,dr,sk03)

        val = val + ejs(3)*rts(3)**6*sk03

        rhoj = -rts(4)
        call sk0(rhoj,dr,sk04)

        val = val + ejs(4)*rts(4)**6*sk04

        rhoj = -rts(5)
        call sk0(rhoj,dr,sk05)

        val = val + ejs(5)*rts(5)**6*sk05

        val = val/2

        return
        end
c
c
c
        subroutine bilapgphiflex(rts,ejs,dr,val)
c
c  returns Delta^2 G_phi(r): biharmonic of the velocity potential Green's function,
c  rts(j)^5 prefactor, scaled by 1/4, plus oneoveralpha/(2*pi*r) free-space term
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(5),ejs(5),val,zt
        complex *16 h0,h1,sk01,sk02,sk03,sk04,sk05
        complex *16 ima,rhoj,oneoveralpha
        real *8 dr,pi
        logical ilow
        data ima /(0,1)/,pi/0.31415926535897932D+01/

        zt = rts(1)*dr

        rhoj = rts(1)

        ione = 1
        call hank103(zt,h0,h1,ione)
        call sk0(rhoj,dr,sk01)

        val = ejs(1)*rts(1)**5*(-sk01+2*ima*h0)

        rhoj = -rts(2)
        call sk0(rhoj,dr,sk02)

        val = val + ejs(2)*rts(2)**5*sk02

        rhoj = -rts(3)
        call sk0(rhoj,dr,sk03)

        val = val + ejs(3)*rts(3)**5*sk03

        rhoj = -rts(4)
        call sk0(rhoj,dr,sk04)

        val = val + ejs(4)*rts(4)**5*sk04

        rhoj = -rts(5)
        call sk0(rhoj,dr,sk05)

        val = val + ejs(5)*rts(5)**5*sk05

        oneoveralpha = ejs(1)*rts(1)**4+ejs(2)*rts(2)**4 
     &           +ejs(3)*rts(3)**4+ejs(4)*rts(4)**4 
     &           +ejs(5)*rts(5)**4

        val = val/4+oneoveralpha/(2*pi*dr)
        
        return
        end
c
c
c
c
c
        subroutine s3dgphiflex(rts,ejs,dr,val)
c
c  returns the 3D free-space part of G_phi at distance dr:
c  residue-weighted K_0 sum without rts prefactor, scaled by 1/8
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        complex *16 rts(5),ejs(5),val,zt
        complex *16 h0,h1,sk01,sk02,sk03,sk04,sk05
        complex *16 ima,rhoj,oneoveralpha
        real *8 dr,pi
        logical ilow
        data ima /(0,1)/,pi/0.31415926535897932D+01/

        zt = rts(1)*dr

        rhoj = rts(1)

        ione = 1
        call hank103(zt,h0,h1,ione)
        call sk0(rhoj,dr,sk01)

        val = ejs(1)*(-sk01+2*ima*h0)

        rhoj = -rts(2)
        call sk0(rhoj,dr,sk02)

        val = val + ejs(2)*sk02

        rhoj = -rts(3)
        call sk0(rhoj,dr,sk03)

        val = val + ejs(3)*sk03

        rhoj = -rts(4)
        call sk0(rhoj,dr,sk04)

        val = val + ejs(4)*sk04

        rhoj = -rts(5)
        call sk0(rhoj,dr,sk05)

        val = val + ejs(5)*sk05

        val = val/8

        return
        end
c
c
c
