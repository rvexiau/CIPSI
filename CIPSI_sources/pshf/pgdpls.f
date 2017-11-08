      subroutine pgdpls
       implicit real*8 (a-h,o-z)
      common/uncp/a(3),b(3),c(3),ga,gb,gc,am,bm,cm,
     1 pt,alfa,alfb,alfc,ctta,r,f1,f2,f3,f4,f5,f6,f7,f8,f9,xbs,bs(10),
     2 xg(64),hg(64),rxg(40),h1(40),h23(40),h456(40),h78(40),h9(40),
     3 fgam(13),rab,l,n,k,inum,nbs
      common/pmat/psdint(225)
      dimension id(6),jd(6)
      data id/1,2,3,1,1,2/,jd/1,2,3,2,3,3/
       do 5 i=1,6
    5 psdint(i)=0.d0
c......
      if(l.lt.0) go to 500
      go to (100,200,300,400),k
c......
  100 call intgab
      do 150 m=1,6
      i=id(m)
      j=jd(m)
      aa=a(i)*a(j)
      am2=am*am
      al=2.d0*alfa*am
      if(l-1) 110,120,130
  110 t1=aa
      t2=-2.d0*aa/am
      if(i.eq.j) t2=t2+1.d0/al
      t3=aa/am2
      go to 140
  120 ab=a(i)*b(j)+a(j)*b(i)
      t1=aa*ctta-ab/(ga*bm)
      t2=(-2.d0*aa*ctta+ab/(ga*bm))/am
      if(i.eq.j) t2=t2+ctta/al
      t3=aa*ctta/am2
      go to 140
  130 pl=1.5d0*ctta*ctta-0.5d0
      ab=a(i)*b(j)+a(j)*b(i)
      bb=b(i)*b(j)
      dab=ga*bm
      ama=am*ga
      w=3.d0*ab*ctta/dab
      t1=aa*(pl+2.d0/ama)+3.d0*bb/(dab*dab)-w
      if(i.eq.j) t1=t1-1.d0/(al*ga)
      t2=(w-2.d0*aa*(pl+1.d0/ama))/am
      if(i.eq.j) t2=t2+pl/al
      t3=aa*pl/am2
  140 continue
      dpls=t1*f1+t2*f2+t3*f3
      psdint(m)=dpls
  150 continue
      go to 1000
c......
  200 if(l-1) 210,1000,220
  210 t1=1.d0/(3.d0*r)
      call intga(gb,5)
      do 215 m=1,3
  215 psdint(m)=t1*f1
      go to 1000
  220 call intga(gb,6)
      fx=-f1/(15.d0*r)
      f1=f1/(5.d0*bm*bm)
      do 225 m=1,6
      i=id(m)
      j=jd(m)
      psdint(m)=b(i)*b(j)*f1
      if(m.le.3) psdint(m)=psdint(m)+fx
  225 continue
      go to 1000
c......
  300 if(l.ne.0) go to 1000
      call intga(ga,7)
      am2=am*am
      al=2.d0*alfa*am
      fx=f2/al
      f2=-f2*2.d0/am
      f3=f3/am2
      fy=f1+f2+f3
      do 310 m=1,6
      i=id(m)
      j=jd(m)
      psdint(m)=a(i)*a(j)*fy
      if(m.le.3) psdint(m)=psdint(m)+fx
  310 continue
      go to 1000
c......
  400 if(l.ne.0) go to 1000
      dpls=fgam(5-n)/(6.0d0*r)
      do 410 m=1,3
  410 psdint(m)=dpls
      go to 1000
c......
  500 go to (600,700,600,900),k
  600 if(cm.lt.1.d-8) go to 650
      s=1.d0/cm
      call intga(gc,18)
      f3=f3*s*s
      f2=f2*s
      fx=f2*s*0.5d0/rab
      do 610 m=1,6
      i=id(m)
      j=jd(m)
      psdint(m)=a(i)*a(j)*f1-(a(i)*c(j)+a(j)*c(i))*f2+c(i)*c(j)*f3
      if(m.le.3) psdint(m)=psdint(m)+fx
  610 continue
      go to 1000
  650 dxpt=dexp(-pt)
      f1=0.5d0*fgam(3-n)*dxpt
      fx=fgam(5-n)*dxpt/(6.0d0*r)
      do 660 m=1,6
      i=id(m)
      j=jd(m)
      psdint(m)=a(i)*a(j)*f1
      if(m.le.3) psdint(m)=psdint(m)+fx
  660 continue
      go to 1000
  700 if(cm.lt.1.d-8) go to 750
      call intga(gc,19)
      s=1.d0/cm
      fx=f1*s*0.5d0/rab
      f2=f2*s*s
      do 710 m=1,6
      i=id(m)
      j=jd(m)
      psdint(m)=c(i)*c(j)*f2
      if(m.le.3) psdint(m)=psdint(m)+fx
  710 continue
      go to 1000
  750 dpls=fgam(5-n)*dexp(-pt)/(6.0d0*r)
      do 760 m=1,3
      psdint(m)=dpls
  760 continue
      go to 1000
c......
  900 dpls=fgam(5-n)/(6.0d0*r)
      do 910 m=1,3
      psdint(m)=dpls
  910 continue
c......
 1000 continue
      return
      end
