      subroutine pgppls
       implicit real*8 (a-h,o-z)
      common/uncp/a(3),b(3),c(3),ga,gb,gc,am,bm,cm,
     1 pt,alfa,alfb,alfc,ctta,r,f1,f2,f3,f4,f5,f6,f7,f8,f9,xbs,bs(10),
     2 xg(64),hg(64),rxg(40),h1(40),h23(40),h456(40),h78(40),h9(40),
     3 fgam(13),rab,l,n,k,inum,nbs
      common /pmat/ psdint(225)
c......
      do 5 i=1,3
    5 psdint(i)=0.d0
c......
      if(l.lt.0) go to 500
      go to (100,200,300,1000),k
c......
  100 call intgab
      do 150 i=1,3
      if(l-1) 110,120,130
  110 t1=-a(i)
      t2=a(i)/am
      go to 140
  120 t1=b(i)/(ga*bm)-a(i)*ctta
      t2=a(i)*ctta/am
      go to 140
  130 pl=1.5d0*ctta*ctta-0.5d0
      aa=a(i)/am
      bb=b(i)/bm
      t1=(3.0d0*bb*ctta-aa)/ga -a(i)*pl
      t2=aa*pl
  140 continue
      ppls=t1*f1+t2*f2
      psdint(i)=ppls
  150 continue
      go to 1000
c......
  200 if(l.ne.1) go to 1000
      call intga(gb,2)
      fx=f1/(3.d0*bm)
      do 210 i=1,3
      ppls=b(i)*fx
      psdint(i)=ppls
  210 continue
      go to 1000
c......
  300 if(l.ne.0) go to 1000
      call intga(ga,3)
      fx=-f1+f2/am
      do 310 i=1,3
      ppls=a(i)*fx
      psdint(i)=ppls
  310 continue
      go to 1000
c......
  500 go to (600,700,600,1000),k
  600 if(cm.lt.1.d-8) go to 650
      call intga(gc,14)
      fx=f2/cm
      do 610 i=1,3
      ppls=-a(i)*f1+c(i)*fx
      psdint(i)=ppls
  610 continue
      go to 1000
  650 fx=-0.5d0*fgam(3-n)*dexp(-pt)
      do 660 i=1,3
  660 psdint(i)=fx*a(i)
      go to 1000
c......
  700 if(cm.lt.1.d-8) go to 1000
      call intga(gc,15)
      fx=f1/cm
      do 710 i=1,3
  710 psdint(i)=fx*c(i)
c......
 1000 continue
      return
      end
