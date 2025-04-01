      implicit real *8 (a-h,o-z)
      complex *16 ima,z0,rhoj0,val0,z1,rhoj1,val1
      data ima /(0,1)/

      z0 = 1.7d0
      rhoj0 = 1.0d0+ima

      z1 = 2.1d0
      rhoj1 = 2.0d0-0.3d0*ima

      call sk0(rhoj0,z0,val0)
      call sk0(rhoj1,z1,val1)

      open(unit=33,file='print_test_sk0.txt')
      write(33,'(2(2x,e22.16))'),val0,val1
      close(33)

      stop
      end
