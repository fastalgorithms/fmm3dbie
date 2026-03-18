

      subroutine el3d_elastlet_mindlin_vec(nd, src, ndt, tar, ndd, &
        dpars, ndz, zpars, ndi, ipars, vals)
!
!  This subroutine computes the greens function tensor 
!  corresponding to the mindlin solution. 
!
!  Input arguments:
!    - nd: integer *8
!        number of kernels returned. Must be atmost 9.
!    - src: real *8(12)
!        source information
!    - ndt: integer *8
!        must be at least 3
!    - tar: real *8(ndt)
!        target information
!    - ndd: integer *8
!        number of real parameters, must at least be 2
!    - dpars: real *8(ndd)
!        dpars(1) - Lame parameter lambda
!        dpars(2) - Lame parameter mu
!    - ndz: integer *8
!        zpars not used
!    - zpars: complex *16(ndz)
!        zpars not used
!    - ndi: integer *8
!        ipars not used
!    - ipars: integer *8(ndi)
!        ipars not used
!
!  Output arguments:
!    - val: real *8(3,3)
!        the Mindlin Greens function 
!
!
      implicit none
      integer *8, intent(in) :: nd, ndt, ndd, ndz, ndi
      real *8, intent(in) :: src(*), tar(ndt)
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer *8, intent(in) :: ipars(ndi)
      real *8, intent(out) :: vals(3,3)

      real *8 dlam, dmu, h
      real *8 rx, ry, rz, dnx, dny, dnz, rfac, bfac, r
      real *8 r3, rdn, tx, ty, tz, rp, rp2, da
      real *8 over4pi
      data over4pi/0.07957747154594767d0/

      dlam = dpars(1)
      dmu = dpars(2)
      da = (dlam + dmu)/(dlam + 2*dmu)
      
      rx = tar(1)-src(1)
      ry = tar(2)-src(2)
      rz = tar(3)-src(3)
        
      dnx = src(10)
      dny = src(11)
      dnz = src(12)

      rfac = over4pi/dmu
      bfac = (1-da)/da*rfac

      r  = sqrt(rx**2+ry**2+rz**2)
      r3 = r**3 
      rdn = dnx*rx+dny*ry+dnz*rz

      tx = rx - rdn*dnx
      ty = ry - rdn*dny
      tz = rz - rdn*dnz
        
      rp = r-rdn
      rp2= rp**2

      vals(1,1) = -tx*tx/rp2/r - (dnx*tx-tx*dnx+rdn*dnx*dnx)/r/rp
      vals(2,1) = -ty*tx/rp2/r - (dny*tx-ty*dnx+rdn*dny*dnx)/r/rp
      vals(3,1) = -tz*tx/rp2/r - (dnz*tx-tz*dnx+rdn*dnz*dnx)/r/rp
      vals(1,2) = -tx*ty/rp2/r - (dnx*ty-tx*dny+rdn*dnx*dny)/r/rp
      vals(2,2) = -ty*ty/rp2/r - (dny*ty-ty*dny+rdn*dny*dny)/r/rp
      vals(3,2) = -tz*ty/rp2/r - (dnz*ty-tz*dny+rdn*dnz*dny)/r/rp
      vals(1,3) = -tx*tz/rp2/r - (dnx*tz-tx*dnz+rdn*dnx*dnz)/r/rp
      vals(2,3) = -ty*tz/rp2/r - (dny*tz-ty*dnz+rdn*dny*dnz)/r/rp
      vals(3,3) = -tz*tz/rp2/r - (dnz*tz-tz*dnz+rdn*dnz*dnz)/r/rp

      vals(1,1) = vals(1,1) + 1/rp
      vals(2,2) = vals(2,2) + 1/rp
      vals(3,3) = vals(3,3) + 1/rp
      vals = bfac*vals

      vals(1,1) = vals(1,1) + rfac*rx*rx/r3
      vals(2,1) = vals(2,1) + rfac*ry*rx/r3
      vals(3,1) = vals(3,1) + rfac*rz*rx/r3
      vals(1,2) = vals(1,2) + rfac*rx*ry/r3
      vals(2,2) = vals(2,2) + rfac*ry*ry/r3
      vals(3,2) = vals(3,2) + rfac*rz*ry/r3
      vals(1,3) = vals(1,3) + rfac*rx*rz/r3
      vals(2,3) = vals(2,3) + rfac*ry*rz/r3
      vals(3,3) = vals(3,3) + rfac*rz*rz/r3

      vals(1,1) = vals(1,1) + rfac/r
      vals(2,2) = vals(2,2) + rfac/r
      vals(3,3) = vals(3,3) + rfac/r


      return
      end
!
!
!
!

      subroutine el3d_elastlet_dn_mindlin_vec(nd, src, ndt, tar, ndd, &
        dpars, ndz, zpars, ndi, ipars, vals)
!
!  This subroutine computes the gradient of the greens function tensor 
!  corresponding to the mindlin solution in the normal source direction 
!
!  Input arguments:
!    - nd: integer *8
!        number of kernels returned. Must be atmost 9.
!    - src: real *8(12)
!        source information
!    - ndt: integer *8
!        must be at least 3
!    - tar: real *8(ndt)
!        target information
!    - ndd: integer *8
!        number of real parameters, must at least be 2
!    - dpars: real *8(ndd)
!        dpars(1) - Lame parameter lambda
!        dpars(2) - Lame parameter mu
!    - ndz: integer *8
!        zpars not used
!    - zpars: complex *16(ndz)
!        zpars not used
!    - ndi: integer *8
!        ipars not used
!    - ipars: integer *8(ndi)
!        ipars not used
!
!  Output arguments:
!    - val: real *8(3,3)
!        the Mindlin Greens function 
!
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: ndt, ndd, ndz, ndi
      real *8, intent(in) :: src(*), tar(ndt)
      real *8, intent(in) :: dpars(ndd)
      real *8, intent(out) :: vals(3,3)

      real *8 dlam, dmu
      real *8 over4pi
      data over4pi/0.07957747154594767d0/

      dlam = dpars(1)
      dmu = dpars(2)
      da = (dlam + dmu)/(dlam + 2*dmu)
      
      rx = tar(1)-src(1)
      ry = tar(2)-src(2)
      rz = tar(3)-src(3)
        
      dnx = src(10)
      dny = src(11)
      dnz = src(12)

      rfac = over4pi/dmu
      bfac = (1-da)/da*rfac

      r  = sqrt(rx**2+ry**2+rz**2)
      r3 = r**3 
      rdn = dnx*rx+dny*ry+dnz*rz

      tx = rx - rdn*dnx
      ty = ry - rdn*dny
      tz = rz - rdn*dnz
      
      rp = r-rdn
      rp2= rp**2

      den = 1/r/rp
      dd1 = (2*r-rdn)/r3/rp2
      dd2 = 1/r3
      dd3 = 1/r/rp

      vals(1,1) = -tx*tx*dd1 - (dnx*tx-tx*dnx+rdn*dnx*dnx)*dd2 - &
                     dnx*dnx*den
      vals(2,1) = -ty*tx*dd1 - (dny*tx-ty*dnx+rdn*dny*dnx)*dd2 - &
                     dny*dnx*den
      vals(3,1) = -tz*tx*dd1 - (dnz*tx-tz*dnx+rdn*dnz*dnx)*dd2 - &
                     dnz*dnx*den
      vals(1,2) = -tx*ty*dd1 - (dnx*ty-tx*dny+rdn*dnx*dny)*dd2 - &
                     dnx*dny*den
      vals(2,2) = -ty*ty*dd1 - (dny*ty-ty*dny+rdn*dny*dny)*dd2 - &
                     dny*dny*den
      vals(3,2) = -tz*ty*dd1 - (dnz*ty-tz*dny+rdn*dnz*dny)*dd2 - &
                     dnz*dny*den
      vals(1,3) = -tx*tz*dd1 - (dnx*tz-tx*dnz+rdn*dnx*dnz)*dd2 - &
                     dnx*dnz*den
      vals(2,3) = -ty*tz*dd1 - (dny*tz-ty*dnz+rdn*dny*dnz)*dd2 - &
                     dny*dnz*den
      vals(3,3) = -tz*tz*dd1 - (dnz*tz-tz*dnz+rdn*dnz*dnz)*dd2 - &
                     dnz*dnz*den

      vals(1,1) = vals(1,1) + dd3
      vals(2,2) = vals(2,2) + dd3
      vals(3,3) = vals(3,3) + dd3
      vals = bfac*vals

      r5 = r**5
      vals(1,1) = vals(1,1) -3*rdn*rfac*rx*rx/r5 + &
                    rfac*(dnx*rx+rx*dnx)/r3 
      vals(2,1) = vals(2,1) -3*rdn*rfac*ry*rx/r5 + &
                    rfac*(dny*rx+ry*dnx)/r3  
      vals(3,1) = vals(3,1) -3*rdn*rfac*rz*rx/r5 + &
                    rfac*(dnz*rx+rz*dnx)/r3
      vals(1,2) = vals(1,2) -3*rdn*rfac*rx*ry/r5 + &
                    rfac*(dnx*ry+rx*dny)/r3
      vals(2,2) = vals(2,2) -3*rdn*rfac*ry*ry/r5 + &
                    rfac*(dny*ry+ry*dny)/r3
      vals(3,2) = vals(3,2) -3*rdn*rfac*rz*ry/r5 + &
                    rfac*(dnz*ry+rz*dny)/r3
      vals(1,3) = vals(1,3) -3*rdn*rfac*rx*rz/r5 + &
                    rfac*(dnx*rz+rx*dnz)/r3
      vals(2,3) = vals(2,3) -3*rdn*rfac*ry*rz/r5 + &
                    rfac*(dny*rz+ry*dnz)/r3
      vals(3,3) = vals(3,3) -3*rdn*rfac*rz*rz/r5 + &
                    rfac*(dnz*rz+rz*dnz)/r3

      vals(1,1) = vals(1,1) -rdn*rfac/r3 
      vals(2,2) = vals(2,2) -rdn*rfac/r3
      vals(3,3) = vals(3,3) -rdn*rfac/r3


      return
      end
!
!
!
!
!
!
      subroutine el3d_elastlet_string_mindlin_vec(nd, src, ndt, tar, &
        ndd, dpars, ndz, zpars, ndi, ipars, vals)
!
!  This subroutine computes the greens function tensor 
!  corresponding to a string of mindlin solution in the normal source direction.
!
!  Let G^{M,3} denote the Green's function with the Mindlin correction, then
!  this subroutine returns G^{M,3}(x,y,n) - G^{M,3}(x,y+hn,n)
!  + hn \cdot \nabla_{y} G^{M,3}(x, y+hn, n), here n is the normal
!  at the source point.
!
!  Input arguments:
!    - nd: integer *8
!        number of kernels returned. Must be atmost 9.
!    - src: real *8(12)
!        source information
!    - ndt: integer *8
!        must be at least 3
!    - tar: real *8(ndt)
!        target information
!    - ndd: integer *8
!        number of real parameters, must at least be 3
!    - dpars: real *8(ndd)
!        dpars(1) - Lame parameter lambda
!        dpars(2) - Lame parameter mu
!        dpars(3) - h (see definition above)
!    - ndz: integer *8
!        zpars not used
!    - zpars: complex *16(ndz)
!        zpars not used
!    - ndi: integer *8
!        ipars not used
!    - ipars: integer *8(ndi)
!        ipars not used
!
!  Output arguments:
!    - val: real *8(3,3)
!        the string of Mindlin Greens function 
!
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: ndt, ndd, ndz, ndi
      real *8, intent(in) :: src(*), tar(ndt)
      real *8, intent(in) :: dpars(ndd)
      real *8, intent(out) :: vals(3,3)

      real *8 dlam, dmu
      real *8 over4pi
      real *8 v1(3,3), v2(3,3), v3(3,3)
      real *8 srctmp(12)
      data over4pi/0.07957747154594767d0/

      dlam = dpars(1)
      dmu = dpars(2)
      da = (dlam + dmu)/(dlam + 2*dmu)

      h = dpars(3)

      srctmp(1:12) = src(1:12)
      srctmp(1) = srctmp(1) + h*src(10)
      srctmp(2) = srctmp(2) + h*src(11)
      srctmp(3) = srctmp(3) + h*src(12)

      call el3d_elastlet_mindlin_vec(nd, src, ndt, tar, ndd, &
        dpars, ndz, zpars, ndi, ipars, v1)

      call el3d_elastlet_mindlin_vec(nd, srctmp, ndt, tar, ndd, &
        dpars, ndz, zpars, ndi, ipars, v2)
        
      call el3d_elastlet_dn_mindlin_vec(nd, srctmp, ndt, tar, ndd, &
        dpars, ndz, zpars, ndi, ipars, v3)

      vals = v1 - v2 - h*v3
      
      return
      end
!
!
!
!
!
      subroutine el3d_elastlet_mindlin_stress_vec(nd, src, ndt, tar, &
        ndd, dpars, ndz, zpars, ndi, ipars, vals)
!
!  This subroutine computes the stress tensor for 
!  greens function tensor corresponding to the mindlin solution. 
!
!  Input arguments:
!    - nd: integer *8
!        number of kernels returned. Must be atmost 27.
!    - src: real *8(12)
!        source information
!    - ndt: integer *8
!        must be at least 3
!    - tar: real *8(ndt)
!        target information
!    - ndd: integer *8
!        number of real parameters, must at least be 2
!    - dpars: real *8(ndd)
!        dpars(1) - Lame parameter lambda
!        dpars(2) - Lame parameter mu
!    - ndz: integer *8
!        zpars not used
!    - zpars: complex *16(ndz)
!        zpars not used
!    - ndi: integer *8
!        ipars not used
!    - ipars: integer *8(ndi)
!        ipars not used
!
!  Output arguments:
!    - val: real *8(3,3,3)
!        stress tensor associated with mindlin Green's function 
!
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: ndt, ndd, ndz, ndi
      real *8, intent(in) :: src(*), tar(ndt)
      real *8, intent(in) :: dpars(ndd)
      real *8, intent(out) :: vals(3,3,3)

      real *8 dlam, dmu
      real *8 rs(3), ns(3), qmat(3,3), eye(3,3), ts(3) 
      real *8 over4pi
      data over4pi/0.07957747154594767d0/

      dlam = dpars(1)
      dmu = dpars(2)
      da = (dlam + dmu)/(dlam + 2*dmu)
      
      rs(1) = tar(1)-src(1)
      rs(2) = tar(2)-src(2)
      rs(3) = tar(3)-src(3)

      ns(1) = src(10)
      ns(2) = src(11)
      ns(3) = src(12)

      r  = sqrt(rs(1)**2+rs(2)**2+rs(3)**2)
      r2 = r**2
      r3 = r**3 
      r5 = r**5 
      rdn = ns(1)*rs(1)+ns(2)*rs(2)+ns(3)*rs(3)

      rfac = over4pi 
      bfac = (1-da)/da*rfac

      ts(1) = rs(1) - rdn*ns(1)
      ts(2) = rs(2) - rdn*ns(2)
      ts(3) = rs(3) - rdn*ns(3)

      qmat(1:3,1:3) = 0

      do i=1,3
        do j=1,3
          qmat(i,j) = -ns(i)*ns(j)
        enddo
      enddo

      qmat(1,1) = 1+qmat(1,1)
      qmat(2,2) = 1+qmat(2,2)
      qmat(3,3) = 1+qmat(3,3)

      eye(1:3,1:3) = 0
      eye(1,1) = 1
      eye(2,2) = 1
      eye(3,3) = 1

      rp = r-rdn
      rp2 = rp**2
      rp3 = rp**3

      do k=1,3
        do j=1,3
          do i=1,3

            rih = ts(i)
            rjh = ts(j)
            dni = ns(i)
            dnj = ns(j)
            dnk = ns(k)
            ri = rs(i)
            rj = rs(j)
            rk = rs(k)
            qij = qmat(i,j)
            qik = qmat(i,k)
            qjk = qmat(j,k)
            eij = eye(i,j)

            t1 = 4*rih*rjh*rk/r2/rp3+2*rih*rjh*rk/r3/rp2
            t2 = -2*qij/r*(1/rp2-1/r2)*rk
            t3 = 2*qij/rp2*dnk
            t4 = -4*rih*rjh*dnk/r/rp3
            t5 = -2*qik*rjh/r/rp2-2*qjk*rih/r/rp2

            s1 = -6*rfac*ri*rj*rk/r5

            vals(i,j,k) = s1 + bfac*(t1 + t2 + t3 + t4 + t5)
          enddo
        enddo
      enddo

      return
      end
!
!
!
!
!

      subroutine el3d_elastlet_dn_mindlin_stress_vec(nd, src, ndt, tar, &
        ndd, dpars, ndz, zpars, ndi, ipars, vals)
!
!  This subroutine computes the stress tensor for 
!  gradient of the greens function tensor corresponding to the mindlin solution
!  in the normal direction.
!
!  Input arguments:
!    - nd: integer *8
!        number of kernels returned. Must be atmost 27.
!    - src: real *8(12)
!        source information
!    - ndt: integer *8
!        must be at least 3
!    - tar: real *8(ndt)
!        target information
!    - ndd: integer *8
!        number of real parameters, must at least be 2
!    - dpars: real *8(ndd)
!        dpars(1) - Lame parameter lambda
!        dpars(2) - Lame parameter mu
!    - ndz: integer *8
!        zpars not used
!    - zpars: complex *16(ndz)
!        zpars not used
!    - ndi: integer *8
!        ipars not used
!    - ipars: integer *8(ndi)
!        ipars not used
!
!  Output arguments:
!    - val: real *8(3,3,3)
!        stress tensor associated with mindlin Green's function 
!
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: ndt, ndd, ndz, ndi
      real *8, intent(in) :: src(*), tar(ndt)
      real *8, intent(in) :: dpars(ndd)
      real *8, intent(out) :: vals(3,3,3)

      real *8 dlam, dmu
      real *8 rs(3), ns(3), qmat(3,3), eye(3,3), ts(3) 
      real *8 over4pi
      data over4pi/0.07957747154594767d0/

      dlam = dpars(1)
      dmu = dpars(2)
      da = (dlam + dmu)/(dlam + 2*dmu)
      

      rs(1) = tar(1)-src(1)
      rs(2) = tar(2)-src(2)
      rs(3) = tar(3)-src(3)

      ns(1) = src(10)
      ns(2) = src(11)
      ns(3) = src(12)

      r  = sqrt(rs(1)**2 + rs(2)**2 + rs(3)**2)
      r2 = r**2
      r3 = r**3 
      r5 = r**5 
      r7 = r**7
      rdn = ns(1)*rs(1) + ns(2)*rs(2) + ns(3)*rs(3)

      rfac = over4pi 
      bfac = (1-da)/da*rfac

      ts(1) = rs(1) - rdn*ns(1)
      ts(2) = rs(2) - rdn*ns(2)
      ts(3) = rs(3) - rdn*ns(3)

      qmat(1:3,1:3) = 0

      do i=1,3
        do j=1,3
          qmat(i,j) = -ns(i)*ns(j)
        enddo
      enddo

      qmat(1,1) = 1+qmat(1,1)
      qmat(2,2) = 1+qmat(2,2)
      qmat(3,3) = 1+qmat(3,3)

      eye(1:3,1:3) = 0
      eye(1,1) = 1
      eye(2,2) = 1
      eye(3,3) = 1

      rp = r-rdn
      rp2= rp**2
      rp3= rp**3
      rp4= rp**4
      rp5= rp**5

      do k=1,3
        do j=1,3
          do i=1,3

            rih = ts(i)
            rjh = ts(j)
            dni = ns(i)
            dnj = ns(j)
            dnk = ns(k)
            ri = rs(i)
            rj = rs(j)
            rk = rs(k)
            qij = qmat(i,j)
            qik = qmat(i,k)
            qjk = qmat(j,k)
            eij = eye(i,j)

            t1 = 4*rih*rjh*rk/r2/rp3+2*rih*rjh*rk/r3/rp2
            t2 = -2*qij/r*(1/rp2-1/r2)*rk
            t3 =  2*qij/rp2*dnk
            t4 = -4*rih*rjh*dnk/r/rp3
            t5 = -2*qik*rjh/r/rp2-2*qjk*rih/r/rp2

            s1 = -6*rfac*(ri*rj*dnk+ri*dnj*rk+dni*rj*rk)/r5 + &
                    30*rfac*rdn*ri*rj*rk/r7
            t14 = -2*rih*rjh*dnk*(rp+2*r)/r3/rp3 + &
                     2*rih*rjh*rk*(2*r2+3*r*rp+3*rp2)/r5/rp3
            t235= 2*(qij*dnk)/r/rp2 - &
                    2*(qik*rjh+qjk*rih+qij*rk)*(rp+r)/r3/rp2 + &
                    2*qij*dnk/r3-6*qij*rk*rdn/r5

            vals(i,j,k) = s1 + bfac*t14 + bfac*t235 

          enddo
        enddo
      enddo

      return
      end
!
!
!
!
      subroutine el3d_elastlet_string_mindlin_normalstress_vec(nd, src, ndt, & 
        tar, ndd, dpars, ndz, zpars, ndi, ipars, vals)
!
!  This subroutine computes the normal stress of the greens function tensor 
!  corresponding to a string of mindlin solution in the normal source direction.
!
!  Let G^{M,3} denote the Green's function with the Mindlin correction, then
!  this subroutine returns nt.\Sigma(G^{M,3}(x,y,n) - G^{M,3}(x,y+hn,n)
!  + hn \cdot \nabla_{y} G^{M,3}(x, y+hn, n)), here n is the normal
!  at the source point, and nt is the target normal
!
!  Input arguments:
!    - nd: integer *8
!        number of kernels returned. Must be atmost 9.
!    - src: real *8(12)
!        source information
!    - ndt: integer *8
!        must be at least 12
!    - tar: real *8(ndt)
!        target information
!    - ndd: integer *8
!        number of real parameters, must at least be 3
!    - dpars: real *8(ndd)
!        dpars(1) - Lame parameter lambda
!        dpars(2) - Lame parameter mu
!        dpars(3) - h (see definition above)
!    - ndz: integer *8
!        zpars not used
!    - zpars: complex *16(ndz)
!        zpars not used
!    - ndi: integer *8
!        ipars not used
!    - ipars: integer *8(ndi)
!        ipars not used
!
!  Output arguments:
!    - val: real *8(3,3)
!        the normal stress of the string of Mindlin Greens function 
!
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: ndt, ndd, ndz, ndi
      real *8, intent(in) :: src(*), tar(ndt)
      real *8, intent(in) :: dpars(ndd)
      real *8, intent(out) :: vals(3,3)

      real *8 dlam, dmu
      real *8 over4pi
      real *8 v1(3,3,3), v2(3,3,3), v3(3,3,3), v4(3,3,3)
      real *8 srctmp(12)
      data over4pi/0.07957747154594767d0/

      dlam = dpars(1)
      dmu = dpars(2)
      da = (dlam + dmu)/(dlam + 2*dmu)

      h = dpars(3)

      srctmp(1:12) = src(1:12)
      srctmp(1) = srctmp(1) + h*src(10)
      srctmp(2) = srctmp(2) + h*src(11)
      srctmp(3) = srctmp(3) + h*src(12)

      call el3d_elastlet_mindlin_stress_vec(nd, src, ndt, tar, ndd, &
        dpars, ndz, zpars, ndi, ipars, v1)

      call el3d_elastlet_mindlin_stress_vec(nd, srctmp, ndt, tar, ndd, &
        dpars, ndz, zpars, ndi, ipars, v2)
        
      call el3d_elastlet_dn_mindlin_stress_vec(nd, srctmp, ndt, tar, ndd, &
        dpars, ndz, zpars, ndi, ipars, v3)
      
      v4 = v1 - v2 - h*v3

      vals(1:3,1:3) = 0
      do jj =1,3
        do ii = 1,3
          do kk = 1,3
            vals(ii,jj) = vals(ii,jj) + v4(kk,ii,jj)*tar(9+kk)
          enddo
        enddo
      enddo

      
      return
      end
!
!
!
!
      subroutine el3d_elastlet_mindlin_normalstress_vec(nd, src, ndt, & 
        tar, ndd, dpars, ndz, zpars, ndi, ipars, vals)
!
!  This subroutine computes the normal stress of the greens function tensor 
!  corresponding to the mindlin solution in the normal source direction.
!
!  Let G^{M,3} denote the Green's function with the Mindlin correction
!
!  Input arguments:
!    - nd: integer *8
!        number of kernels returned. Must be atmost 9.
!    - src: real *8(12)
!        source information
!    - ndt: integer *8
!        must be at least 12
!    - tar: real *8(ndt)
!        target information
!    - ndd: integer *8
!        number of real parameters, must at least be 3
!    - dpars: real *8(ndd)
!        dpars(1) - Lame parameter lambda
!        dpars(2) - Lame parameter mu
!        dpars(3) - h (see definition above)
!    - ndz: integer *8
!        zpars not used
!    - zpars: complex *16(ndz)
!        zpars not used
!    - ndi: integer *8
!        ipars not used
!    - ipars: integer *8(ndi)
!        ipars not used
!
!  Output arguments:
!    - val: real *8(3,3)
!        the normal stress of the string of Mindlin Greens function 
!
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: ndt, ndd, ndz, ndi
      real *8, intent(in) :: src(*), tar(ndt)
      real *8, intent(in) :: dpars(ndd)
      real *8, intent(out) :: vals(3,3)

      real *8 dlam, dmu
      real *8 over4pi
      real *8 v1(3,3,3), v2(3,3,3), v3(3,3,3), v4(3,3,3)
      real *8 srctmp(12)
      data over4pi/0.07957747154594767d0/

      dlam = dpars(1)
      dmu = dpars(2)
      da = (dlam + dmu)/(dlam + 2*dmu)

      h = dpars(3)

      srctmp(1:12) = src(1:12)

      call el3d_elastlet_mindlin_stress_vec(nd, src, ndt, tar, ndd, &
        dpars, ndz, zpars, ndi, ipars, v1)

      vals(1:3,1:3) = 0
      do jj =1,3
        do ii = 1,3
          do kk = 1,3
            vals(ii,jj) = vals(ii,jj) + v1(kk,ii,jj)*tar(9+kk)
          enddo
        enddo
      enddo

      
      return
      end
!
!
!
!
