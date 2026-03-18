!
!
!  This file contains the parametrizations of several
!  analytic charts
!
!
! 
      subroutine wtorus_eval(s, t, radii, scales, nosc, xyz, dxyzdst)
!
!  This subrotine returns the chart info for a wiggly torus
!  whose parametrization is given by
!  
!  x(s,t) = a1*r(s,t)*cos(s)
!  y(s,t) = a2*r(s,t)*sin(s)
!  z(s,t) = a3*r1*sin(t)
!
!  r(s,t) = r2 + r1*cos(t) + r3*cos(nosc*s)
!
!  radii = (r1, r2, r3)
!  scales = (a1, a2, a3)
!
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      real *8 s, t, radii(3), scales(3)
      integer *8 nosc
      real *8 xyz(3), dxyzdst(3,2)
      
      real *8 rr, drrds, drrdt

      rr = radii(2) + radii(1)*cos(t) + radii(3)*cos(nosc*s)
      xyz(1) = scales(1)*rr*cos(s)
      xyz(2) = scales(2)*rr*sin(s)
      xyz(3) = scales(3)*radii(1)*sin(t)

      drrds = -nosc*radii(3)*sin(nosc*s)
      drrdt = -radii(1)*sin(t)

      dxyzdst(1,1) = scales(1)*(drrds*cos(s) - rr*sin(s))
      dxyzdst(2,1) = scales(2)*(drrds*sin(s) + rr*cos(s))
      dxyzdst(3,1) = 0

      dxyzdst(1,2) = scales(1)*drrdt*cos(s)
      dxyzdst(2,2) = scales(2)*drrdt*sin(s)
      dxyzdst(3,2) = scales(3)*radii(1)*cos(t)
      
      return
      end



      subroutine xyz_tensor_fourier_eval(s, t, coefs, m, nfp, scales, &
        xyz, dxyzdst)
!
!  This subroutine is the chart for a double fourier toroidal
!  discretization
!
!  \hat(x) = \sum_{i=1}^{2m+1} \sum_{j=1} x_{ij} b_{i}(s) b_{j}(nfp*t)
!  \hat(y) = \sum_{i=1}^{2m+1} \sum_{j=1} y_{ij} b_{i}(s) b_{j}(nfp*t)
!  \hat(z) = \sum_{i=1}^{2m+1} \sum_{j=1} z_{ij} b_{i}(s) b_{j}(nfp*t)
!
!  x(s,t) = (\hat(x) \cos(t) - \hat(y) \sin(t))*scales(1)
!  y(s,t) = (\hat(x) \sin(t) + \hat(y) \cos(t))*scales(2)
!  z(s,t) = \hat(z)*scales(3)

      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      real *8 :: xyz(3), dxyzdst(3,2), coefs(2*m+1,2*m+1,3)
      real *8 :: hatxyz(3), dhatxyzds(3), dhatxyzdt(3)
      real *8 :: bis(2*m+1), bjs(2*m+1), bdis(2*m+1), bdjs(2*m+1)
      real *8 :: scales(3)
      complex *16 zmuls, zmult, ima, zfacs, zfact
      integer *8 nfp
      data ima/(0.0d0,1.0d0)/

      xyz(1) = 0
      xyz(2) = 0
      xyz(3) = 0

      hatxyz(1:3) = 0

      dxyzdst(1,1) = 0
      dxyzdst(2,1) = 0
      dxyzdst(3,1) = 0

      dhatxyzds(1:3) = 0

      dxyzdst(1,2) = 0
      dxyzdst(2,2) = 0
      dxyzdst(3,2) = 0

      dhatxyzdt(1:3) = 0
  
      ct = cos(t)
      st = sin(t)
      
      c4t = cos(nfp*t)
      s4t = sin(nfp*t)
      zmuls = cos(s) + ima*sin(s)
      zmult = c4t + ima*s4t


      zfacs = zmuls
      zfact = zmult
  
      bis(1) = 1
      bdis(1) = 0
  
      bjs(1) = 1
      bdjs(1) = 0

      do i=1,m
        bis(i+1) = real(zmuls)
        bis(i+1+m) = imag(zmuls)

        bdis(i+1) = -i*imag(zmuls)
        bdis(i+m+1) = i*real(zmuls)
     
        bjs(i+1) = real(zmult)
        bjs(i+1+m) = imag(zmult)

        bdjs(i+1) = -nfp*i*imag(zmult)
        bdjs(i+1+m) = nfp*i*real(zmult)


        zmuls = zmuls*zfacs
        zmult = zmult*zfact
      enddo
!
!
      do j = 1,2*m+1
        do i = 1,2*m+1
          hatxyz(1) = hatxyz(1) + coefs(i,j,1)*bis(i)*bjs(j)
          hatxyz(2) = hatxyz(2) + coefs(i,j,2)*bis(i)*bjs(j)
          hatxyz(3) = hatxyz(3) + coefs(i,j,3)*bis(i)*bjs(j)

          dhatxyzds(1) = dhatxyzds(1) + coefs(i,j,1)*bdis(i)*bjs(j)
          dhatxyzds(2) = dhatxyzds(2) + coefs(i,j,2)*bdis(i)*bjs(j)
          dhatxyzds(3) = dhatxyzds(3) + coefs(i,j,3)*bdis(i)*bjs(j)

          dhatxyzdt(1) = dhatxyzdt(1) + coefs(i,j,1)*bis(i)*bdjs(j)
          dhatxyzdt(2) = dhatxyzdt(2) + coefs(i,j,2)*bis(i)*bdjs(j)
          dhatxyzdt(3) = dhatxyzdt(3) + coefs(i,j,3)*bis(i)*bdjs(j)
        enddo
      enddo

      xyz(1) = (hatxyz(1)*ct - hatxyz(2)*st)*scales(1)
      xyz(2) = (hatxyz(1)*st + hatxyz(2)*ct)*scales(2)
      xyz(3) = hatxyz(3)*scales(3)
  
      dxyzdst(1,1) = (dhatxyzds(1)*ct - dhatxyzds(2)*st)*scales(1) 
      dxyzdst(2,1) = (dhatxyzds(1)*st + dhatxyzds(2)*ct)*scales(2) 
      dxyzdst(3,1) = dhatxyzds(3)*scales(3)

      dxyzdst(1,2) = -(hatxyz(1)*st + hatxyz(2)*ct)
      dxyzdst(1,2) = dxyzdst(1,2) + (dhatxyzdt(1)*ct - dhatxyzdt(2)*st)
      dxyzdst(1,2) = dxyzdst(1,2)*scales(1)

      dxyzdst(2,2) = (hatxyz(1)*ct - hatxyz(2)*st)
      dxyzdst(2,2) = dxyzdst(2,2) + (dhatxyzdt(1)*st + dhatxyzdt(2)*ct)
      dxyzdst(2,2) = dxyzdst(2,2)*scales(2)

      dxyzdst(3,2) = dhatxyzdt(3)*scales(3)

      return
      end



      subroutine torus_curve_tube_eval(s,t,dpars,xyz,dxyzdst)
!
!  This subroutine parametrizes a tubular surface around a torus curve.
!
!  The center curve is
!
!      gamma(s) = ((rmaj + rmin*cos(polfq*s))*cos(torfq*s),
!                  (rmaj + rmin*cos(polfq*s))*sin(torfq*s),
!                   rmin*sin(polfq*s))
!
!  which lies on a torus with major radius rmaj and minor radius rmin.
!  The parameters polfq and torfq control how the curve winds in the
!  poloidal and toroidal directions, respectively.
!
!  Note: polfq/torfq must be rational for a closed center curve.
!
!  The surface is obtained by attaching a circular cross-section of
!  radius rcross along the curve:
!
!      x(s,t) = (rho(s) + rcross*cos(polfq*s)*cos(t))*cos(torfq*s)
!               + rcross*a1(s)*sin(t)
!
!      y(s,t) = (rho(s) + rcross*cos(polfq*s)*cos(t))*sin(torfq*s)
!               + rcross*b1(s)*sin(t)
!
!      z(s,t) =  rmin*sin(polfq*s) + rcross*sin(polfq*s)*cos(t)
!               - rcross*g1(s)*sin(t)
!
!  where
!
!      rho(s) = rmaj + rmin*cos(polfq*s)
!
!      a1(s) = ( torfq*rho(s)*sin(polfq*s)*cos(torfq*s)
!                - rmin*polfq*sin(torfq*s) ) / l(s)
!
!      b1(s) = ( rmin*polfq*cos(torfq*s)
!                + torfq*rho(s)*sin(polfq*s)*sin(torfq*s) ) / l(s)
!
!      g1(s) = ( torfq*rho(s)*cos(polfq*s) ) / l(s)
!
!      l(s)  = sqrt(rmin^2*polfq^2 + torfq^2*rho(s)^2)
!
!  The parameter s runs along the torus curve and t is the angular
!  coordinate around the tubular cross-section.
!
!  input:
!    s,t       - surface parameters
!    dpars(*)  - shape parameters, stored as
!                dpars(1) = rmaj
!                dpars(2) = rmin
!                dpars(3) = rcross      
!                dpars(4) = polfq
!                dpars(5) = torfq
!
!  output:
!    xyz(1:3)      - position vector
!    dxyzdst(:,1)  - derivative with respect to s
!    dxyzdst(:,2)  - derivative with respect to t
!

      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)

      real *8 :: s,t
      real *8 :: xyz(3), dxyzdst(3,2), dpars(*)

      real *8 :: rmaj,rmin,polfq,torfq,rcross
      real *8 :: rho,rhos,l,ls
      real *8 :: cc,sc,cd,sd,ct,st
      real *8 :: na,nb,ng
      real *8 :: nas,nbs,ngs
      real *8 :: aa,bb,gg
      real *8 :: aas,bbs,ggs

      rmaj   = dpars(1)
      rmin   = dpars(2)
      rcross = dpars(3)
      polfq  = dpars(4)
      torfq  = dpars(5)


      cc = cos(polfq*s)
      sc = sin(polfq*s)
      cd = cos(torfq*s)
      sd = sin(torfq*s)
      ct = cos(t)
      st = sin(t)

      rho  = rmaj + rmin*cc
      rhos = -rmin*polfq*sc

      l  = sqrt(rmin*rmin*polfq*polfq + torfq*torfq*rho*rho)
      ls = torfq*torfq*rho*rhos/l

      na = torfq*rho*sc*cd - rmin*polfq*sd
      nb = rmin*polfq*cd + torfq*rho*sc*sd
      ng = torfq*rho*cc

      aa = na/l
      bb = nb/l
      gg = ng/l

      nas = torfq*(rhos*sc + rho*polfq*cc)*cd &
           - torfq*torfq*rho*sc*sd &
           - rmin*polfq*torfq*cd

      nbs = -rmin*polfq*torfq*sd &
           + torfq*(rhos*sc + rho*polfq*cc)*sd &
           + torfq*torfq*rho*sc*cd

      ngs = torfq*(rhos*cc - rho*polfq*sc)

      aas = (nas*l - na*ls)/(l*l)
      bbs = (nbs*l - nb*ls)/(l*l)
      ggs = (ngs*l - ng*ls)/(l*l)

      xyz(1) = (rho + rcross*cc*ct)*cd + rcross*aa*st
      xyz(2) = (rho + rcross*cc*ct)*sd + rcross*bb*st
      xyz(3) = rmin*sc + rcross*sc*ct - rcross*gg*st

      dxyzdst(1,1) = (rhos - rcross*polfq*sc*ct)*cd &
                   - torfq*(rho + rcross*cc*ct)*sd &
                   + rcross*aas*st

      dxyzdst(2,1) = (rhos - rcross*polfq*sc*ct)*sd &
                   + torfq*(rho + rcross*cc*ct)*cd &
                   + rcross*bbs*st

      dxyzdst(3,1) = rmin*polfq*cc &
                   + rcross*polfq*cc*ct &
                   - rcross*ggs*st

      dxyzdst(1,2) = -rcross*cc*cd*st + rcross*aa*ct
      dxyzdst(2,2) = -rcross*cc*sd*st + rcross*bb*ct
      dxyzdst(3,2) = -rcross*sc*st - rcross*gg*ct

      return
      
      end


      subroutine fourier_curve_tube_eval(s,t,dpars,xyz,dxyzdst)
!
!  This subroutine parametrizes the surface
!
!      x(s,t) = r*sin(3*s) / (alpha + cos(t))
!
!      y(s,t) = r*(sin(s) + 2*sin(2*s))
!     &         / (alpha + cos(t + phi))
!
!      z(s,t) = (r/8)*(cos(s) - 2*cos(2*s))
!     &         * (c0 + cos(t))
!     &         * (alpha + cos(t + phi))
!
!  The original example corresponds to
!
!      r     = 5
!      alpha = 1.9
!      phi   = 2*pi/3
!      c0    = 2
!
!  input:
!    s,t       - surface parameters
!    dpars(*)  - shape parameters, stored as
!                dpars(1) = r
!                dpars(2) = alpha
!                dpars(3) = phi
!                dpars(4) = c0
!
!  output:
!    xyz(1:3)      - position vector
!    dxyzdst(:,1)  - derivative with respect to s
!    dxyzdst(:,2)  - derivative with respect to t
!

      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)

      real *8 :: s,t
      real *8 :: xyz(3), dxyzdst(3,2), dpars(*)

      real *8 :: r,alpha,phi,c0
      real *8 :: d1,d2,tp
      real *8 :: s1,s2,s3,c1,c2,c3
      real *8 :: ct,st,ctp,stp
      real *8 :: q,qs

      r     = dpars(1)
      alpha = dpars(2)
      phi   = dpars(3)
      c0    = dpars(4)

      s1  = sin(s)
      s2  = sin(2.0d0*s)
      s3  = sin(3.0d0*s)
      c1  = cos(s)
      c2  = cos(2.0d0*s)
      c3  = cos(3.0d0*s)

      ct  = cos(t)
      st  = sin(t)

      tp  = t + phi
      ctp = cos(tp)
      stp = sin(tp)

      d1 = alpha + ct
      d2 = alpha + ctp

      q  = c1 - 2.0d0*c2
      qs = -s1 + 4.0d0*s2

      xyz(1) = r*s3/d1
      xyz(2) = r*(s1 + 2.0d0*s2)/d2
      xyz(3) = (r/8.0d0)*q*(c0 + ct)*d2

      dxyzdst(1,1) = 3.0d0*r*c3/d1
      dxyzdst(2,1) = r*(c1 + 4.0d0*c2)/d2
      dxyzdst(3,1) = (r/8.0d0)*qs*(c0 + ct)*d2

      dxyzdst(1,2) = r*s3*st/(d1*d1)
      dxyzdst(2,2) = r*(s1 + 2.0d0*s2)*stp/(d2*d2)
      dxyzdst(3,2) = -(r/8.0d0)*q*(st*d2 + (c0 + ct)*stp)

      return
      end

      subroutine cruller_eval(s,t,radii,scales,xyz,dxyzdst)
!
!  This subroutine returns the chart info for a cruller-like torus
!  whose parametrization is given by
!
!    h(s,t) = rmean + ramp*cos(freqs*s + freqt*t)
!
!    x(s,t) = sx*(rmaj + h(s,t)*cos(t))*cos(s)
!    y(s,t) = sy*(rmaj + h(s,t)*cos(t))*sin(s)
!    z(s,t) = sz*h(s,t)*sin(t)
!
!  radii  = (rmaj, rmean, ramp, freqs, freqt)
!  scales = (sx, sy, sz)
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)

      real *8 s, t, radii(5), scales(3)
      real *8 xyz(3), dxyzdst(3,2)

      real *8 rmaj, rmean, ramp, freqs, freqt
      real *8 h, dhds, dhdt
      real *8 ct, st, cs, ss
      real *8 rr, drrds, drrdt

      rmaj  = radii(1)
      rmean = radii(2)
      ramp  = radii(3)
      freqs = radii(4)
      freqt = radii(5)

      cs = cos(s)
      ss = sin(s)
      ct = cos(t)
      st = sin(t)

      h = rmean + ramp*cos(freqs*s + freqt*t)

      rr = rmaj + h*ct

      xyz(1) = scales(1)*rr*cs
      xyz(2) = scales(2)*rr*ss
      xyz(3) = scales(3)*h*st

      dhds = -ramp*freqs*sin(freqs*s + freqt*t)
      dhdt = -ramp*freqt*sin(freqs*s + freqt*t)

      drrds = dhds*ct
      drrdt = dhdt*ct - h*st

      dxyzdst(1,1) = scales(1)*(drrds*cs - rr*ss)
      dxyzdst(2,1) = scales(2)*(drrds*ss + rr*cs)
      dxyzdst(3,1) = scales(3)*dhds*st

      dxyzdst(1,2) = scales(1)*drrdt*cs
      dxyzdst(2,2) = scales(2)*drrdt*ss
      dxyzdst(3,2) = scales(3)*(dhdt*st + h*ct)

      return
      end
