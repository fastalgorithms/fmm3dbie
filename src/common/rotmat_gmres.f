
      subroutine rotmat_gmres(a,b,c,s)
      implicit real *8 (a-h,o-z)
      real *8 a,b,c,s

      if(b.eq.0) then
        c = 1
        s = 0
      else if(abs(b).gt.abs(a)) then
        temp = a/b
        s = 1.0d0/sqrt(1.0d0+temp**2)
        c = temp*s
      else
        temp = b/a
        c = 1.0d0/sqrt(1.0d0+temp**2)
        s = temp*c
      endif

      return
      end

          



      subroutine zrotmat_gmres(a,b,c,s)
      implicit real *8 (a-h,o-z)
      complex *16 a,b,c,s,temp

      if(b.eq.0) then
        c = 1
        s = 0
      else if(abs(b).gt.abs(a)) then
        temp = a/b
        s = 1.0d0/sqrt(1.0d0+temp**2)
        c = temp*s
      else
        temp = b/a
        c = 1.0d0/sqrt(1.0d0+temp**2)
        s = temp*c
      endif

      return
      end

          


