      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      ntests = 1
      call test_find_near(i1)

      ntests = ntests + 6
      call test_adap_quad_self(i2)

      ntests = ntests + 2
      call test_mixed_patch_quad(i3)

      ntests = ntests + 3
      call test_near_quad_sub(i4)

      ntests = ntests + 2
      call test_near_corr_bel(i5)

      print *, i1, i2, i3, i4, i5

      nsuccess = i1 + i2 + i3 + i4 + i5

      open(unit=33,file='../../print_testres.txt',access='append')
      write(33,'(a,i2,a,i2,a)') 'Successfully completed ',nsuccess,
     1  ' out of ',ntests,' in quadratures testing suite'
      close(33)
      
      
      

      stop
      end
