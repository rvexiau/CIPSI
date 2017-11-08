      subroutine pgpplp
       implicit real*8 (a-h,o-z)
      logical*1 iandj
      common/uncp/a(3),b(3),c(3),ga,gb,gc,am,bm,cm,
     1 pt,alfa,alfb,alfc,ctta,r,f1,f2,f3,f4,f5,f6,f7,f8,f9,xbs,bs(10),
     2 xg(64),hg(64),rxg(40),h1(40),h23(40),h456(40),h78(40),h9(40),
     3 fgam(13),rab,l,n,k,inum,nbs
      common /pmat/ psdint(225),pspt(225),iandj
c......
      do 5 i=1,9
    5 psdint(i)=0.d0
      maxj=3
      if(l.lt.0) go to 500
      go to (100,200,300,400),k
c......
  100 call intgab
      ij=0
      do 150 i=1,3
      if(iandj)maxj=i
      do 150 j=1,maxj
      ij=ij+1
      aa=a(i)*a(j)
      bb=b(i)*b(j)
      ab=a(i)*b(j)
      if(l-1) 110,120,130
  110 t1=ab
      t2=-ab/am
      t3=-ab/bm
      t4=ab/(am*bm)
      go to 140
  120 w=ab*ctta
      da=aa/(am*gb)
      db=bb/(bm*ga)
      t1=w-da-db
      if(i.eq.j) t1=t1+1.d0/(ga*gb*r)
      t2=(da-w)/am
      t3=(db-w)/bm
      t4=w/(am*bm)
      go to 140
  130 pl=1.5d0*ctta*ctta-0.5d0
      ba=b(i)*a(j)
      da=aa/(am*gb)
      db=bb/(bm*ga)
      ca=ab/am
      cb=ab/bm
      ama=1.d0/(am*ga)
      bmb=1.d0/(bm*gb)
      t1=ab*(pl+ama+bmb)+(3.d0*ba-2.d0*ab)*ama*bmb-3.d0*(da+db)*ctta
      if(i.eq.j) t1=t1+3.d0*ctta/(r*ga*gb)
      t2=3.d0*da*ctta/am-ca*(pl+bmb)
      t3=3.d0*db*ctta/bm-cb*(pl+ama)
      t4=ca*pl/bm
  140 continue
      pplp=t1*f1+t2*f2+t3*f3+t4*f4
      psdint(ij)=pplp
  150 continue
      go to 1000
c......
  200 if(l.ne.1) go to 1000
      call intga(gb,4)
      f1=f1*0.3333333333333333d0
      f2=f2*0.3333333333333333d0/(bm*bm)
      fx=f1/(gb*r)
      f1=-f1/bm
      fy=f1+f2
      ij=0
      do 210 i=1,3
      if(iandj) maxj=i
      do 210 j=1,maxj
      ij=ij+1
      pplp=b(i)*b(j)*fy
      if(i.eq.j) pplp=pplp+fx
      psdint(ij)=pplp
  210 continue
      go to 1000
c......
  300 if(l.ne.1) go to 1000
      call intga(ga,4)
      f1=f1*0.3333333333333333d0
      f2=f2*0.3333333333333333d0/(am*am)
      fx=f1/(ga*r)
      f1=-f1/am
      fy=f1+f2
      ij=0
      do 310 i=1,3
      if(iandj) maxj=i
      do 310 j=1,maxj
      ij=ij+1
      pplp=a(i)*a(j)*fy
      if(i.eq.j) pplp=pplp+fx
      psdint(ij)=pplp
  310 continue
      go to 1000
c......
  400 if(l.ne.1) go to 1000
      pplp=fgam(5-n)/(18.0d0*r)
      ij=0
      do 410 i=1,3
      if(iandj) maxj=i
      do 410 j=1,maxj
      ij=ij+1
      if(i.ne.j) go to 410
      psdint(ij)=pplp
  410 continue
      go to 1000
c......
  500 go to (600,700,800,900),k
  600 if(cm.lt.1.d-8) go to 650
      call intga(gc,16)
      s=1.d0/cm
      f3=f3*s*s
      f2=f2*s
      fx=f2*0.5d0/rab
      ij=0
      do 610 i=1,3
      if(iandj) maxj=i
      do 610 j=1,maxj
      ij=ij+1
      psdint(ij)=a(i)*b(j)*f1-(a(i)*c(j)+c(i)*b(j))*f2+c(i)*c(j)*f3
      if(i.eq.j) psdint(ij)=psdint(ij)+fx
  610 continue
      go to 1000
  650 dxpt=dexp(-pt)
      f1=0.5d0*fgam(3-n)*dxpt
      fx=fgam(5-n)*dxpt/(6.d0*r)
      ij=0
      do 660 i=1,3
      if(iandj) maxj=i
      do 660 j=1,maxj
      ij=ij+1
      psdint(ij)=a(i)*b(j)*f1
      if(i.eq.j) psdint(ij)=psdint(ij)+fx
  660 continue
      go to 1000
c......
  700 if(cm.lt.1.d-8) go to 750
      call intga(gc,17)
      s=1.d0/cm
      f1=f1*s
      f2=f2*s*s
      fx=f1*0.5d0/rab
      ij=0
      do 710 i=1,3
      if(iandj) maxj=i
      do 710 j=1,maxj
      ij=ij+1
      psdint(ij)=-c(i)*b(j)*f1+c(i)*c(j)*f2
      if(i.eq.j) psdint(ij)=psdint(ij)+fx
  710 continue
      go to 1000
  750 pplp=fgam(5-n)*dexp(-pt)/(6.d0*r)
      ij=0
      do 760 i=1,3
      if(iandj) maxj=i
      do 760 j=1,maxj
      ij=ij+1
      if(i.ne.j) go to 760
      psdint(ij)=pplp
  760 continue
      go to 1000
c......
  800 if(cm.lt.1.d-8) go to 750
      call intga(gc,17)
      s=1.d0/cm
      f1=f1*s
      f2=f2*s*s
      fx=f1*0.5d0/rab
      ij=0
      do 810 i=1,3
      if(iandj) maxj=i
      do 810 j=1,maxj
      ij=ij+1
      psdint(ij)=-a(i)*c(j)*f1+c(i)*c(j)*f2
      if(i.eq.j) psdint(ij)=psdint(ij)+fx
  810 continue
      go to 1000
c......
  900 pplp=fgam(5-n)/(6.d0*r)
      ij=0
      do 910 i=1,3
      if(iandj) maxj=i
      do 910 j=1,maxj
      ij=ij+1
      if(i.ne.j) go to 910
      psdint(ij)=pplp
  910 continue
c......
 1000 continue
      return
      end
