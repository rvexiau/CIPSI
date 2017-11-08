      subroutine dipint
c
c     ----- gauss-hermite quadrature using minimum point formula -----
c
      implicit real*8 (a-h,o-z)
      common/intdip/xint0,yint0,zint0,xintx,yinty,zintz,t,x0,y0,z0,
     1 xi,yi,zi,xj,yj,zj,ni,nj
      common/hermit/h1,h2(2),h3(3),h4(4),h5(5),h6(6)
      common/wermit/w1,w2(2),w3(3),w4(4),w5(5),w6(6)
      dimension h(21),w(21),min(6),max(6)
      equivalence (h(1),h1),(w(1),w1)
      data min /1,2,4,7,11,16/
      data max /1,3,6,10,15,21/
      data zero /0.0d+00/
      xint0=zero
      yint0=zero
      zint0=zero
      xintx=zero
      yinty=zero
      zintz=zero
      npts=(ni+nj+1)/2+1
      imin=min(npts)
      imax=max(npts)
      do 13 i=imin,imax
      dum=w(i)
      px=dum
      py=dum
      pz=dum
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
    1 px=px*ax
      py=py*ay
      pz=pz*az
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
   12 continue
      xint0=xint0+px
      yint0=yint0+py
      zint0=zint0+pz
      xintx=xintx+px*ptx
      yinty=yinty+py*pty
      zintz=zintz+pz*ptz
   13 continue
      return
      end
