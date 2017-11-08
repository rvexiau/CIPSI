      subroutine denint
      implicit real*8 (a-h,o-z)
      common/intden/xint,yint,zint,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      data one /1.0d+00/
      xint=one
      yint=one
      zint=one
      ptxi=x0-xi
      ptyi=y0-yi
      ptzi=z0-zi
      ptxj=x0-xj
      ptyj=y0-yj
      ptzj=z0-zj
      go to (5,4,3,2,1),ni
    1 xint=xint*ptxi
      yint=yint*ptyi
      zint=zint*ptzi
    2 xint=xint*ptxi
      yint=yint*ptyi
      zint=zint*ptzi
    3 xint=xint*ptxi
      yint=yint*ptyi
      zint=zint*ptzi
    4 xint=xint*ptxi
      yint=yint*ptyi
      zint=zint*ptzi
    5 go to (10,9,8,7,6),nj
    6 xint=xint*ptxj
      yint=yint*ptyj
      zint=zint*ptzj
    7 xint=xint*ptxj
      yint=yint*ptyj
      zint=zint*ptzj
    8 xint=xint*ptxj
      yint=yint*ptyj
      zint=zint*ptzj
    9 xint=xint*ptxj
      yint=yint*ptyj
      zint=zint*ptzj
   10 continue
      return
      end
