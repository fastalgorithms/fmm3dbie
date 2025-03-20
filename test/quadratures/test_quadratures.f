      implicit real *8 (a-h,o-z)
      ntests = 1
      call test_find_near(i1)

      ntests = ntests + 6
      call test_adap_quad_self(i2)

      print *, i1, i2

      nsuccess = i1 + i2

      open(unit=33,file='../../print_testres.txt',access='append')
      write(33,'(a,i2,a,i2,a)') 'Successfully completed ',nsuccess,
     1  ' out of ',ntests,' in quadratures testing suite'
      close(33)
      
      
      

      stop
      end
