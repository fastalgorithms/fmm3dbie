c (C) The Simons Foundation, Inc., Michael O'Neil, and Felipe Vico 
c
c
c  This file contains the following user-callable routines
c  
c    - rotmat_gmres: 
c        compute the sin and cos such that tan = b/a, with stability
c        close to 0,1
c    - zrotmat_gmres:
c        complex version of rotmat_gmres
c



      subroutine rotmat_gmres(a,b,c,s)
c
c-----------------------
c  Given a,b, compute sin(theta), cos(theta),
c  such that tan(theta) = b/a, note, the routine
c  implements a stabler version of computing
c  b/sqrt(a^2+b^2) and a/sqrt(b^2+a^2)
c
c  Input arguments:
c  
c    - a: double precision
c        cos scaling
c    - b: double precision
c        sin scaling
c
c  Output arguments:
c    - c: double precision
c        cos(theta)
c    - s: double precision
c        sin(theta)
c        
c
c-----------------------
c
      implicit real *8 (a-h,o-z)
      real *8 a,b,c,s

      if(a.eq.0) then
        c = 0
        s = 1
      else
        temp = b/a
        c = 1.0d0/sqrt(1.0d0+temp**2)
        s = temp*c

      endif

      return
      end

          



      subroutine zrotmat_gmres(a,b,c,s)
c-----------------------
c  Given a,b, compute sin(theta), cos(theta),
c  such that tan(theta) = b/a, note, the routine
c  implements a stabler version of computing
c  b/sqrt(a^2+b^2) and a/sqrt(b^2+a^2).
c
c  This routine is the complex version of rotmat_gmres
c
c  Input arguments:
c  
c    - a: double complex
c        cos scaling
c    - b: double complex 
c        sin scaling
c
c  Output arguments:
c    - c: double complex 
c        cos(theta)
c    - s: double complex 
c        sin(theta)
c        
c
c-----------------------
      implicit real *8 (a-h,o-z)
      complex *16 a,b,c,s,temp

      if(a.eq.0) then
        s = 1
        c = 0
      else
        temp = b/a
        c = 1.0d0/sqrt(1.0d0+abs(temp)**2)
        s = temp*c
      endif

      return
      end

          


