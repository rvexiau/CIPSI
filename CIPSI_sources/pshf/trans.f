      subroutine trans(nn)
      implicit real*8 (a-h,o-z)
c
c      calculate the coordinates (xnew,ynew,znew) of the transform
c      of the point (xold,yold,zold) under the transformation t
c
      common/transf/xold,yold,zold,xnew,ynew,znew,xp,yp,zp
      common/symmat/t(432)
      xnew=xold*t(nn+1)+yold*t(nn+2)+zold*t(nn+3)
      ynew=xold*t(nn+4)+yold*t(nn+5)+zold*t(nn+6)
      znew=xold*t(nn+7)+yold*t(nn+8)+zold*t(nn+9)
      return
      end
