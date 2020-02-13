c
c (c) Zydrunas Gimbutas 
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging routines and the beginning of
c       the geometry processing routines 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine dot_prod3d(x,y,d)
        implicit real *8 (a-h,o-z)
        dimension x(3),y(3)
c
c       d = x \dot y
c
        d=x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
c
        return
        end
c
c
c
c
c
        subroutine cross_prod3d(x,y,z)
        implicit real *8 (a-h,o-z)
        dimension x(3),y(3),z(3)
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
        implicit real *8 (a-h,o-z)
        complex *16 x(3),y(3),z(3)
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
        subroutine dot_cross_prod3d(x,y,z,d)
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
        implicit real *8 (a-h,o-z)
        dimension x(3),y(3),z(3),w(3),t(3)
c
c       w = x \cross (y \cross z)
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
        implicit real *8 (a-h,o-z)
        dimension x(3),y(3),z(3),w(3),t(3)
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
