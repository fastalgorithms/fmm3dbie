      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      character *1 ttype

      ifp = 0
      ttype = 'F'
      istot = 0
      ntot = 0
      nquadmax = 5000
      do istrat = 1,3
        do intype=1,2
          do ipoly = 0,1

            isuccess = 0
            call dquadintrouts_testing(istrat,ifp,intype,ipoly,ttype,
     1         nquadmax,isuccess)
            istot = istot + isuccess
            isuccess = 0
            call quadintrouts_testing(istrat,ifp,intype,ipoly,ttype,
     1         nquadmax,isuccess)
            istot = istot + isuccess
            ntot = ntot + 2 
          enddo
        enddo
      enddo
c
c  test that the adaptive integration routines correctly
c  reallocate the tree arrays when nquadmax is too small
c  to begin with
c
      istrat = 1
      nquadmax = 10
      do intype=1,2
        do ipoly = 0,1

          isuccess = 0
          call dquadintrouts_testing(istrat,ifp,intype,ipoly,ttype,
     1       nquadmax,isuccess)
          istot = istot + isuccess
          isuccess = 0
          call quadintrouts_testing(istrat,ifp,intype,ipoly,ttype,
     1       nquadmax,isuccess)
          istot = istot + isuccess
          ntot = ntot + 2
        enddo
      enddo

      open(unit=33,file='../../print_testres.txt',access='append')
      write(33,'(a,i2,a,i2,a)') 'Successfully completed ',istot,
     1  ' out of ',ntot,' in quad_routs testing suite'
      close(33)
      
      return
      end
          

