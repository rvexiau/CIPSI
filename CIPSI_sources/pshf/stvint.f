      subroutine stvint
      implicit real*8 (a-h,o-z)
      common/stv/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      common/hermit/h1,h2(2),h3(3),h4(4),h5(5),h(6)
      common/wermit/w1,w2(2),w3(3),w4(4),w5(5),w(6)
      xint=0.0d+00
      yint=0.0d+00
      zint=0.0d+00
      do 13 i=1,6
      px=1.0d+00
      py=1.0d+00
      pz=1.0d+00
      dum=h(i)/t
      ptx=dum+x0
      pty=dum+y0
      ptz=dum+z0
      ax=ptx-xi
      ay=pty-yi
      az=ptz-zi
      bx=ptx-xj
      by=pty-yj
      bz=ptz-zj
      go to (5,4,3,2,1),ni
    1 px=   ax
      py=   ay
      pz=   az
    2 px=px*ax
      py=py*ay
      pz=pz*az
    3 px=px*ax
      py=py*ay
      pz=pz*az
    4 px=px*ax
      py=py*ay
      pz=pz*az
    5 go to (12,11,10,9,8,7,6),nj
    6 px=px*bx
      py=py*by
      pz=pz*bz
    7 px=px*bx
      py=py*by
      pz=pz*bz
    8 px=px*bx
      py=py*by
      pz=pz*bz
    9 px=px*bx
      py=py*by
      pz=pz*bz
   10 px=px*bx
      py=py*by
      pz=pz*bz
   11 px=px*bx
      py=py*by
      pz=pz*bz
   12 dum=w(i)
      xint=xint+dum*px
      yint=yint+dum*py
      zint=zint+dum*pz
   13 continue
      return
      end
