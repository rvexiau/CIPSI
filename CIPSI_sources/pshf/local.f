      subroutine local(x,y,z,xs,ys,zs)
      implicit real*8 (a-h,o-z)
c
c     calculate the coordinates (xs,ys,zs) of a point in the local
c     frame given the coordinates (x,y,z) in the master frame
c
      common/frame/u1,u2,u3,v1,v2,v3,w1,w2,w3,x0,y0,z0
      xs=u1*(x-x0)+u2*(y-y0)+u3*(z-z0)
      ys=v1*(x-x0)+v2*(y-y0)+v3*(z-z0)
      zs=w1*(x-x0)+w2*(y-y0)+w3*(z-z0)
      return
      end
