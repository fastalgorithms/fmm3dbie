      implicit real *8 (a-h,o-z)
      implicit integer(8) (i-n)


      ntests = 7
      call test_triarouts(i1)

      ntests = ntests + 6
      call test_triarouts_vec(i2)
      nsuccess = i1+i2

      ntests = ntests+ 7
      call test_dtriarouts(i1)

      ntests = ntests + 6
      call test_dtriarouts_vec(i2)

      nsuccess = nsuccess + i1+i2

      ntests = ntests + 1
      call test_koornf_ders(i3)
      nsuccess = nsuccess+i3

      print *,i1,i2,i3

      open(unit=33,file='../../print_testres.txt',access='append')
      write(33,'(a,i2,a,i2,a)') 'Successfully completed ',nsuccess,
     1  ' out of ',ntests,' in tria_routs testing suite'
      close(33)
      
      stop
      end
