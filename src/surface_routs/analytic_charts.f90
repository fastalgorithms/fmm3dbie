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
      real *8 s, t, radii(3), scales(3)
      integer nosc
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



      subroutine xyz_tensor_fourier_eval(s, t, coefs, m, scales, &
        xyz, dxyzdst)
!
!  This subroutine is the chart for a double fourier toroidal
!  discretization
!
!  \hat(x) = \sum_{i=1}^{2m+1} \sum_{j=1} x_{ij} b_{i}(s) b_{j}(t)
!  \hat(y) = \sum_{i=1}^{2m+1} \sum_{j=1} y_{ij} b_{i}(s) b_{j}(t)
!  \hat(z) = \sum_{i=1}^{2m+1} \sum_{j=1} z_{ij} b_{i}(s) b_{j}(t)
!
!  x(s,t) = (\hat(x) \cos(t) - \hat(y) \sin(t))*scales(1)
!  y(s,t) = (\hat(x) \sin(t) + \hat(y) \cos(t))*scales(2)
!  z(s,t) = \hat(z)*scales(3)

      implicit real *8 (a-h,o-z)
      real *8 :: xyz(3), dxyzdst(3,2), coefs(2*m+1,2*m+1,3)
      real *8 :: hatxyz(3), dhatxyzds(3), dhatxyzdt(3)
      real *8 :: bis(2*m+1), bjs(2*m+1), bdis(2*m+1), bdjs(2*m+1)
      real *8 :: scales(3)
      complex *16 zmuls, zmult, ima, zfacs, zfact
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
      zmuls = cos(s) + ima*sin(s)
      zmult = ct + ima*st


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

        bdjs(i+1) = -i*imag(zmult)
        bdjs(i+1+m) = i*real(zmult)


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




