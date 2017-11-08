      subroutine rot
      implicit real*8 (a-h,o-z)
c
c     calculate the coordinates (xp,yp,zp) of a point in the master
c     frame given the coordinates (xnew,ynew,znew) in the local frame
c
      common/transf/xold,yold,zold,xnew,ynew,znew,xp,yp,zp
      common/frame/u1,u2,u3,v1,v2,v3,w1,w2,w3,x0,y0,z0
      xp=x0+u1*xnew+v1*ynew+w1*znew
      yp=y0+u2*xnew+v2*ynew+w2*znew
      zp=z0+u3*xnew+v3*ynew+w3*znew
      return
      end
