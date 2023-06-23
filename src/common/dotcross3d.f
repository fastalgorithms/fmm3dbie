c (c) Zydrunas Gimbutas 
c 
c  This file contains the following user callable routines:
c  
c  * dot_prod3d:
c      compute the dot product of two vectors
c  * zdot_prod3d:
c      computes the dot product of two complex vectors (no conjugation)
c  * dzdot_prod3d:
c      computes the dot product of a real vector and a complex vector
c  * cross_prod3d: 
c      compute the cross product of two vectors
c  * zcross_prod3d:
c      compute the cross product of two complex vectors
c  * dzcross_prod3d:
c      compute the cross product of a real vector and a complex vector 
c  * dot_cross_prod3d:
c      compute x \cdot (y \times z)
c  * cross_cross_prod3d:
c      compute x \times (y \times z)
c  * dot_cross_cross_prod3d
c      compute x \cdot (y \times (z \times w))
c
C  The following routines were written by Felipe Vico and added by Mike O'Neil
C
C  * orthonormalize
c
c
        subroutine dot_prod3d(x,y,d)
c
c------------------
c  This subroutine computes the dot product of two vectors
c
c  Input arguments:
c
c    - x: double precision(3)
c        input vector 1
c    - y: double precision(3)
c        input vector 2
c  
c  Output arguments:
c
c    - d: double precision
c        x \cdot y
c
c----------------
        implicit real *8 (a-h,o-z)
        double precision, intent(in) :: x(3),y(3)
        double precision, intent(out) :: d

        d=x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
c
        return
        end
c
c
c
c
c
        subroutine zdot_prod3d(x,y,d)
c
c------------------
c  This subroutine computes the dot product of two vectors
c
c  Input arguments:
c
c    - x: double complex(3)
c        input vector 1
c    - y: double complex(3)
c        input vector 2
c  
c  Output arguments:
c
c    - d: double complex
c        x \cdot y
c
c----------------
        implicit real *8 (a-h,o-z)
        double complex, intent(in) :: x(3),y(3)
        double complex, intent(out) :: d

        d=x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
c
        return
        end
c
c
c
c
        subroutine dzdot_prod3d(x,y,d)
c
c------------------
c  This subroutine computes the dot product of 
c  a real vector with a complex vector
c
c  Input arguments:
c
c    - x: double precision(3)
c        input vector 1
c    - y: double complex(3)
c        input vector 2
c  
c  Output arguments:
c
c    - d: double complex
c        x \cdot y
c
c----------------
        implicit real *8 (a-h,o-z)
        double precision, intent(in) :: x(3)
        double complex, intent(in) :: y(3)
        double complex, intent(out) :: d

        d=x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
c
        return
        end
c
c
c
c
c
c
        subroutine cross_prod3d(x,y,z)
c------------------
c  This subroutine computes the cross product of two vectors
c
c  Input arguments:
c
c    - x: double precision(3)
c        input vector 1
c    - y: double precision(3)
c        input vector 2
c  
c  Output arguments:
c
c    - z: double precision(3)
c        x \times y
c---------------------
        implicit real *8 (a-h,o-z)
        double precision, intent(in) :: x(3),y(3)
        double precision, intent(out) :: z(3)
c
c       z = x \cross y
c
        z(1)=x(2)*y(3)-x(3)*y(2)
        z(2)=x(3)*y(1)-x(1)*y(3)
        z(3)=x(1)*y(2)-x(2)*y(1)
c
        return
        end
c
c
c
c
c
        subroutine zcross_prod3d(x,y,z)
c------------------
c  This subroutine computes the cross product of two complex vectors
c
c  Input arguments:
c
c    - x: double complex(3)
c        input vector 1
c    - y: double complex(3)
c        input vector 2
c  
c  Output arguments:
c
c    - z: double complex(3)
c        x \times y
c---------------------
        implicit real *8 (a-h,o-z)
        complex *16, intent(in) :: x(3),y(3)
        complex *16, intent(out) :: z(3)
c
c       z = x \cross y
c
        z(1)=x(2)*y(3)-x(3)*y(2)
        z(2)=x(3)*y(1)-x(1)*y(3)
        z(3)=x(1)*y(2)-x(2)*y(1)
c
        return
        end
c
c
c
c
c
        subroutine dzcross_prod3d(x,y,z)
c------------------
c  This subroutine computes the cross product of a real vector
c  with a complex vector
c
c  Input arguments:
c
c    - x: double precision(3)
c        input vector 1
c    - y: double complex(3)
c        input vector 2
c  
c  Output arguments:
c
c    - z: double complex(3)
c        x \times y
c---------------------
        implicit real *8 (a-h,o-z)
        real *8, intent(in) :: x(3)
        complex *16, intent(in) :: y(3)
        complex *16, intent(out) :: z(3)
c
c       z = x \cross y
c
        z(1)=x(2)*y(3)-x(3)*y(2)
        z(2)=x(3)*y(1)-x(1)*y(3)
        z(3)=x(1)*y(2)-x(2)*y(1)
c
        return
        end
c
        subroutine dot_cross_prod3d(x,y,z,d)
c------------------
c  This subroutine computes the scalar triple product 
c
c  Input arguments:
c
c    - x: double precision(3)
c        input vector 1
c    - y: double precision(3)
c        input vector 2
c    - z: double precision(3)
c        input vector 3
c  
c  Output arguments:
c
c    - d: double precision
c        x \cdot (y \times z)
c---------------------
        implicit real *8 (a-h,o-z)
        dimension x(3),y(3),z(3),t(3)
c
c       d = x \dot (y \cross z)
c
        call cross_prod3d(y,z,t)
        call dot_prod3d(x,t,d)
c
        return
        end
c
c
c
c
c
        subroutine cross_cross_prod3d(x,y,z,w)
c------------------
c  This subroutine computes the vector triple product 
c
c  Input arguments:
c
c    - x: double precision(3)
c        input vector 1
c    - y: double precision(3)
c        input vector 2
c    - z: double precision(3)
c        input vector 3
c  
c  Output arguments:
c
c    - w: double precision(3)
c        x \times (y \times z)
c---------------------
        implicit real *8 (a-h,o-z)
        double precision, intent(in) :: x(3),y(3),z(3)
        double precision, intent(out) :: w(3)
        double precision t(3)
c
        call cross_prod3d(y,z,t)
        call cross_prod3d(x,t,w)
c
        return
        end
c
c
c
c
c
        subroutine dot_cross_cross_prod3d(x,y,z,w,d)
c------------------
c  This subroutine computes the dot product with vector triple product 
c
c  Input arguments:
c
c    - x: double precision(3)
c        input vector 1
c    - y: double precision(3)
c        input vector 2
c    - z: double precision(3)
c        input vector 3
c    - w: double precision(3)
c        input vector 4
c  
c  Output arguments:
c
c    - d: double precision(3)
c        x \cdot (y \times (z \times w))
c---------------------
        implicit real *8 (a-h,o-z)
        double precision, intent(in) :: x(3),y(3),z(3),w(3)
        double precision, intent(out) :: d
        double precision t(3)
c
c       d = x \dot (y \cross (z \cross w))
c
        call cross_cross_prod3d(y,z,w,t)
        call dot_prod3d(x,t,d)
c
        return
        end
c
c
c
c
c

        subroutine orthonormalize_all(du,normal,ru,rv,ns)
c!f2py intent(in) du, normal, ns
c!f2py intent(out) ru,rv
        implicit none
c       !List of calling arguments
        integer, intent(in) :: ns
        double precision, intent(in) :: du(3,ns), normal(3,ns)
        double precision, intent(out) :: ru(3,ns), rv(3,ns)

c       !List of local variables
        double precision :: aux
        integer :: j

        do j=1,ns
           call orthonormalize(du(:,j), normal(:,j), ru(:,j), rv(:,j))
         enddo

        return
        end






        subroutine orthonormalize(du, normal, ru, rv)
        implicit none
c
c       This routine computes:
C           ru = du / ||du||
C           rv = normal X ru
c
c       !List of calling arguments
        double precision, intent(in) :: du(3), normal(3)
        double precision, intent(out) :: ru(3), rv(3)
c
c       !List of local variables
        double precision :: aux
c
        aux=sqrt(du(1)**2+du(2)**2+du(3)**2)
        ru(1)=du(1)/aux
        ru(2)=du(2)/aux
        ru(3)=du(3)/aux
c
        call cross_prod3d(normal, ru, rv)
c        call my_cross_v2(normal, ru, rv)
        return
        end






csubroutine my_cross_v2(a, b, c)
cimplicit none
c
c    real ( kind = 8 ), intent(in) :: a(3),b(3)
c    real ( kind = 8 ), intent(out) :: c(3)
c
c        c(1) = a(2) * b(3) - a(3) * b(2)
c        c(2) = a(3) * b(1) - a(1) * b(3)
c        c(3) = a(1) * b(2) - a(2) * b(1)
c
cend subroutine my_cross_v2
