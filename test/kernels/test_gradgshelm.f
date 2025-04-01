      implicit real *8 (a-h,o-z)
      complex *16 ima,rhoj0,rhoj1
      complex *16 gradx, grady
      complex *16 rts(3),ejs(3)
      real *8 dx0,dx1,dy0,dy1
      data ima /(0,1)/

      rts(1) = 2.0d0
      rts(2) = -1.0d0 + ima
      rts(3) = -1.0d0 - ima

      ejs(1) = 0.2d0
      ejs(2) = -0.1d0 + 0.3d0*ima
      ejs(3) = -0.1d0 - 0.3d0*ima

      dx0 = 1.1d0
      dy0 = 2.3d0

      call gradgshelm(rts,ejs,dx0,dy0,gradx,grady)

      open(unit=33,file='print_test_gradgshelm.txt')
      write(33,'(4(2x,e22.16))'),gradx,grady
      close(33)

      stop
      end
