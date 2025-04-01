      implicit real *8 (a-h,o-z)
      complex *16 rts(3),ejs(3)
      complex *16 ima,z0,val
      data ima /(0,1)/

      rts(1) = 2.0d0 
      rts(2) = -1.0d0 + ima
      rts(3) = -1.0d0 - ima

      ejs(1) = 0.2d0
      ejs(2) = -0.1d0 + 0.3d0*ima
      ejs(3) = -0.1d0 - 0.3d0*ima

      z0 = 1.7d0

      call lapgphihelm(rts,ejs,z0,val)

      open(unit=33,file='print_test_lapgphihelm.txt')
      write(33,'(1(2x,e22.16))'),val
      close(33)

      stop
      end
