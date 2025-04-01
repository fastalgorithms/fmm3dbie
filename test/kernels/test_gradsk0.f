      implicit real *8 (a-h,o-z)
      complex *16 ima,rhoj0,rhoj1
      complex *16 sk0x0,sk0y0,sk0x1,sk0y1
      real *8 dx0,dx1,dy0,dy1
      data ima /(0,1)/

      dx0 = 1.1d0
      dy0 = 2.3d0
      rhoj0 = 1.0d0+ima

      dx1 = 2.1d0
      dy1 = 0.7d0
      rhoj1 = 2.0d0-0.3d0*ima

      call gradsk0(rhoj0,dx0,dy0,sk0x0,sk0y0)
      call gradsk0(rhoj1,dx1,dy1,sk0x1,sk0y1)

      open(unit=33,file='print_test_gradsk0.txt')
      write(33,'(4(2x,e22.16))'),sk0x0,sk0y0,sk0x1,sk0y1
      close(33)

      stop
      end
