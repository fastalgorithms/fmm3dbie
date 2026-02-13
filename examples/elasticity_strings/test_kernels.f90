      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      real *8 src(12), tar(12), val_grid(3,3,-2:2,-2:2,-2:2)
      real *8 src_grid(12, -2:2,-2:2,-2:2), uval(3,-2:2,-2:2,-2:2)
      real *8 str(3), uxx(3), uyy(3), uzz(3), uxy(3), uxz(3), uyx(3)
      real *8 uyz(3), uzx(3), uzy(3), dpars(3), ux(3), uy(3), uz(3)
      real *8 ucomp(3), utest(3), utest2(3)
      real *8 amat(3,3)
      complex *16 zpars
      integer *8 ipars, ndd, ndz, ndi
      
      done = 1.0d0
      pi = atan(one)*4
      call prini(6,13)
      
      src(1:12) = 0 
      src(1) = hkrand(0)*2 - 1.0d0
      src(2) = hkrand(0)*2 - 1.0d0
      src(3) = hkrand(0)*2 - 1.0d0
      
      src(10) = hkrand(0)*2 - 1.0d0
      src(11) = hkrand(0)*2 - 1.0d0
      src(12) = hkrand(0)*2 - 1.0d0
      
      rr = sqrt(src(10)**2 + src(11)**2 + src(12)**2)
      src(10:12) = src(10:12)/rr

      dmu = 1.3d0
      dla = 1.4d0
      h = 0.3d0
      dpars(1) = dla
      dpars(2) = dmu
      dpars(3) = h
      da = (dla + dmu)/(dla + 2*dmu)

      tar(1:12) = 0
      tar(1) = hkrand(0)*2 - 1.0d0
      tar(2) = hkrand(0)*2 - 1.0d0
      tar(3) = hkrand(0)*2 - 1.0d0
      
      tar(10) = hkrand(0)*2 - 1.0d0
      tar(11) = hkrand(0)*2 - 1.0d0
      tar(12) = hkrand(0)*2 - 1.0d0

!      tar(10) = 0
!      tar(11) = 0
!      tar(12) = 1


      
      rr = sqrt(tar(10)**2 + tar(11)**2 + tar(12)**2)
      tar(10:12) = tar(10:12)/rr

      str(1) = hkrand(0)
      str(2) = hkrand(0)
      str(3) = hkrand(0)

      str(1:3) = 0
      str(2) = 1

      hh = 1.0d-3

      nd = 9
      ndt = 12
      ndd = 3
      ndi = 0
      ndz = 0

      do ilat = -2,2
        do jlat = -2,2
          do klat = -2,2
            src_grid(1:12,ilat,jlat,klat) = src(1:12)
            src_grid(1,ilat,jlat,klat) = src(1) + ilat*hh
            src_grid(2,ilat,jlat,klat) = src(2) + jlat*hh
            src_grid(3,ilat,jlat,klat) = src(3) + klat*hh
            val_grid(1:3,1:3,ilat,jlat,klat) = 0
            call el3d_elastlet_string_mindlin_vec(nd, &
              src_grid(1,ilat,jlat,klat), ndt, tar, ndd, dpars, ndz, &
              zpars, ndi, ipars, val_grid(1,1,ilat,jlat,klat))
              
              uval(1:3,ilat,jlat,klat) = & 
                  val_grid(1:3,1,ilat,jlat,klat)*str(1) + &      
                  val_grid(1:3,2,ilat,jlat,klat)*str(2) + &      
                  val_grid(1:3,3,ilat,jlat,klat)*str(3)        

          enddo
        enddo
      enddo
      uxy(1:3) = (uval(1:3,1,1,0) - uval(1:3,1,-1,0) - & 
                    uval(1:3,-1,1,0) + uval(1:3,-1,-1,0))/4/hh/hh
      uxz(1:3) = (uval(1:3,1,0,1) - uval(1:3,1,0,-1) - &
                    uval(1:3,-1,0,1) + uval(1:3,-1,0,-1))/4/hh/hh
      uyz(1:3) = (uval(1:3,0,1,1) - uval(1:3,0,1,-1) - &
                    uval(1:3,0,-1,1) + uval(1:3,0,-1,-1))/4/hh/hh
      uyx(1:3) = uxy(1:3)
      uzx(1:3) = uxz(1:3)
      uzy(1:3) = uyz(1:3)

      uxx(1:3) = (uval(1:3,1,0,0) - 2*uval(1:3,0,0,0) + &
            uval(1:3,-1,0,0))/hh/hh
      uyy(1:3) = (uval(1:3,0,1,0) - 2*uval(1:3,0,0,0) + &
            uval(1:3,0,-1,0))/hh/hh
      uzz(1:3) = (uval(1:3,0,0,1) - 2*uval(1:3,0,0,0) + &
            uval(1:3,0,0,-1))/hh/hh
      call prin2('uxx=*', uxx, 3)
      call prin2('uyy=*', uyy, 3)
      call prin2('uzz=*', uzz, 3)
      call prin2('uxy=*', uxy, 3)
      call prin2('uxz=*', uxz, 3)
      call prin2('uyz=*', uyz, 3)
      c1 = (uxx(1) + uxy(2) + uxz(3))*(dmu + dla) + &
            dmu*(uxx(1) + uyy(1) + uzz(1))
      c2 = (uyx(1) + uyy(2) + uyz(3))*(dmu + dla) + &
            dmu*(uxx(2) + uyy(2) + uzz(2))
      c3 = (uzx(1) + uzy(2) + uzz(3))*(dmu + dla) + &
            dmu*(uxx(3) + uyy(3) + uzz(3))
      call prin2('c1=*', c1, 1)
      call prin2('c2=*', c2, 1)
      call prin2('c3=*', c3, 1)

      ux(1:3) = (uval(1:3,1,0,0) - uval(1:3,-1,0,0))/2/hh;
      uy(1:3) = (uval(1:3,0,1,0) - uval(1:3,0,-1,0))/2/hh;
      uz(1:3) = (uval(1:3,0,0,1) - uval(1:3,0,0,-1))/2/hh;

      divu = ux(1) + uy(2) + uz(3)
      amat(1:3,1:3) = 0
      amat(1,1) = divu*dla
      amat(2,2) = divu*dla
      amat(3,3) = divu*dla

      amat(1,1) = amat(1,1) + dmu*(ux(1) + ux(1))
      amat(1,2) = amat(1,2) + dmu*(ux(2) + uy(1))
      amat(1,3) = amat(1,3) + dmu*(ux(3) + uz(1))
      
      amat(2,1) = amat(2,1) + dmu*(uy(1) + ux(2))
      amat(2,2) = amat(2,2) + dmu*(uy(2) + uy(2))
      amat(2,3) = amat(2,3) + dmu*(uy(3) + uz(2))
      
      amat(3,1) = amat(3,1) + dmu*(uz(1) + ux(3))
      amat(3,2) = amat(3,2) + dmu*(uz(2) + uy(3))
      amat(3,3) = amat(3,3) + dmu*(uz(3) + uz(3))

      utest(1:3) = tar(10)*amat(1,1:3) + tar(11)*amat(2,1:3) + &
         tar(12)*amat(3,1:3)
      utest2(1:3) = tar(10)*amat(1:3,1) + tar(11)*amat(1:3,2) + &
         tar(12)*amat(1:3,3)
      call el3d_elastlet_string_mindlin_normalstress_vec(nd, &
              src, ndt, tar, ndd, dpars, ndz, &
              zpars, ndi, ipars, ucomp)
      call prin2('utest=*',utest,3)
      call prin2('utest2=*',utest2,3)
      call prin2('ucomp=*',ucomp,3)
      


      stop
      end
