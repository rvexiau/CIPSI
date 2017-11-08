      subroutine xyzint
      implicit real*8 (a-h,o-z)
      logical*1 n0,n1,m0,m1
      common/setint/i(9),k(9),nimax,njmax,nkmax,nlmax,nmax,mmax
     1,bp01,b00,b10,xcp00,xc00,ycp00,yc00,zcp00,zc00,f00
     2,dxij,dyij,dzij,dxkl,dykl,dzkl
      common/xyz/xint(5625),yint(5625),zint(5625)
      n0=nmax.eq.0
      n1=nmax.le.1
      m0=mmax.eq.0
      m1=mmax.le.1
c
c     ----- i(0,0) -----
c
      i1=i(1)
      xint(i1)=1.0d+00
      yint(i1)=1.0d+00
      zint(i1)=f00
      if(n0.and.m0) return
      if(n0) go to 5
c
c     ----- i(1,0) -----
c
      i2=i(2)
      xint(i2)=xc00
      yint(i2)=yc00
      zint(i2)=zc00*f00
    5 if(m0) go to 10
c
c     ----- i(0,1) -----
c
      k2=k(2)
      i3=i1+k2
      xint(i3)=xcp00
      yint(i3)=ycp00
      zint(i3)=zcp00*f00
      if(n0) go to 10
c
c     ----- i(1,1) -----
c
      i3=i2+k2
      cp10=b00
      xint(i3)=xcp00*xint(i2)+cp10
      yint(i3)=ycp00*yint(i2)+cp10
      zint(i3)=zcp00*zint(i2)+cp10*f00
   10 if(n1) go to 40
      c10=0.0d+00
      i3=i1
      i4=i2
      do 30 n=2,nmax
      c10=c10+b10
c
c     ----- i(n,0) -----
c
      i5=i(n+1)
      xint(i5)=c10*xint(i3)+xc00*xint(i4)
      yint(i5)=c10*yint(i3)+yc00*yint(i4)
      zint(i5)=c10*zint(i3)+zc00*zint(i4)
      if(m0) go to 20
      cp10=cp10+b00
c
c     ----- i(n,1) -----
c
      i3=i5+k2
      xint(i3)=xcp00*xint(i5)+cp10*xint(i4)
      yint(i3)=ycp00*yint(i5)+cp10*yint(i4)
      zint(i3)=zcp00*zint(i5)+cp10*zint(i4)
   20 i3=i4
   30 i4=i5
   40 if(m1) go to 70
      cp01=0.0d+00
      c01=b00
      i3=i1
      i4=i1+k2
      do 60 m=2,mmax
      cp01=cp01+bp01
c
c     ----- i(0,m) -----
c
      i5=i1+k(m+1)
      xint(i5)=cp01*xint(i3)+xcp00*xint(i4)
      yint(i5)=cp01*yint(i3)+ycp00*yint(i4)
      zint(i5)=cp01*zint(i3)+zcp00*zint(i4)
      if(n0) go to 50
      c01=c01+b00
c
c     ----- i(1,m) -----
c
      i3=i2+k(m+1)
      xint(i3)=xc00*xint(i5)+c01*xint(i4)
      yint(i3)=yc00*yint(i5)+c01*yint(i4)
      zint(i3)=zc00*zint(i5)+c01*zint(i4)
   50 i3=i4
   60 i4=i5
   70 if(n1.or.m1) go to 100
c
c     ----- i(n,m) -----
      c01=b00
      k3=k2
      do 90 m=2,mmax
      k4=k(m+1)
      c01=c01+b00
      i3=i1
      i4=i2
      c10=b10
      do 80 n=2,nmax
      i5=i(n+1)
      xint(i5+k4)=c10*xint(i3+k4)+xc00*xint(i4+k4)+c01*xint(i4+k3)
      yint(i5+k4)=c10*yint(i3+k4)+yc00*yint(i4+k4)+c01*yint(i4+k3)
      zint(i5+k4)=c10*zint(i3+k4)+zc00*zint(i4+k4)+c01*zint(i4+k3)
      c10=c10+b10
      i3=i4
   80 i4=i5
   90 k3=k4
  100 if(njmax.eq.0) go to 200
c
c     ----- i(ni,nj,m) -----
c
      m=0
      i5=i(nmax+1)
  110 min=nimax
      km=k(m+1)
  120 n=nmax
      i3=i5+km
  130 i4=i(n)+km
      xint(i3)=xint(i3)+dxij*xint(i4)
      yint(i3)=yint(i3)+dyij*yint(i4)
      zint(i3)=zint(i3)+dzij*zint(i4)
      i3=i4
      n=n-1
      if(n.gt.min) go to 130
      min=min+1
      if(min.lt.nmax) go to 120
      if(nimax.eq.0) go to 160
      i3=25+km+i1
      do 150 nj=1,njmax
      i4=i3
      do 140 ni=1,nimax
      xint(i4)=xint(i4+100)+dxij*xint(i4-25)
      yint(i4)=yint(i4+100)+dyij*yint(i4-25)
      zint(i4)=zint(i4+100)+dzij*zint(i4-25)
  140 i4=i4+125
  150 i3=i3+25
  160 m=m+1
      if(m.le.mmax) go to 110
  200 if(nlmax.eq.0) go to 300
c
c     ----- i(ni,nj,nk,nl) -----
c
      i5=k(mmax+1)
      ia=i1
      ni=0
  210 nj=0
      ib=ia
  220 min=nkmax
  230 m=mmax
      i3=ib+i5
  240 i4=ib+k(m)
      xint(i3)=xint(i3)+dxkl*xint(i4)
      yint(i3)=yint(i3)+dykl*yint(i4)
      zint(i3)=zint(i3)+dzkl*zint(i4)
      i3=i4
      m=m-1
      if(m.gt.min) go to 240
      min=min+1
      if(min.lt.mmax) go to 230
      if(nkmax.eq.0) go to 270
      i3=ib+1
      do 260 nl=1,nlmax
      i4=i3
      do 250 nk=1,nkmax
      xint(i4)=xint(i4+4)+dxkl*xint(i4-1)
      yint(i4)=yint(i4+4)+dykl*yint(i4-1)
      zint(i4)=zint(i4+4)+dzkl*zint(i4-1)
  250 i4=i4+5
  260 i3=i3+1
  270 nj=nj+1
      ib=ib+25
      if(nj.le.njmax) go to 220
      ni=ni+1
      ia=ia+125
      if(ni.le.nimax) go to 210
  300 continue
      return
      end
