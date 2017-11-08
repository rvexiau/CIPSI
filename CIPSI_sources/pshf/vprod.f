      subroutine vprod(vp,x,y)
c
c     vp=x cross y
c
      implicit real*8 (a-h,o-z)
      dimension vp(3),x(3),y(3)
      vp(1)=x(2)*y(3)-x(3)*y(2)
      vp(2)=x(3)*y(1)-x(1)*y(3)
      vp(3)=x(1)*y(2)-x(2)*y(1)
      return
      end
