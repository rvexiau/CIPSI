      subroutine pgspls
       implicit real*8 (a-h,o-z)
      common/uncp/a(3),b(3),c(3),ga,gb,gc,am,bm,cm,
     1 pt,alfa,alfb,alfc,ctta,r,f1,f2,f3,f4,f5,f6,f7,f8,f9,xbs,bs(10),
     2 xg(64),hg(64),rxg(40),h1(40),h23(40),h456(40),h78(40),h9(40),
     3 fgam(13),rab,l,n,k,inum,nbs
      common/pmat/psdint(225)
      equivalence (psdint(1),spls)
c......
      spls=0.d0
      if(l.lt.0) go to 500
      go to (100,200,300,400),k
c......
  100 call intgab
      if(l-1) 110,120,130
  110 spls=f1
      go to 1000
  120 spls=f1*ctta
      go to 1000
  130 spls=(1.5d0*ctta*ctta-0.5d0)*f1
      go to 1000
c......
  200 if(l.ne.0) go to 1000
      call intga(gb,1)
      spls=f1
      go to 1000
c......
  300 if(l.ne.0) go to 1000
      call intga(ga,1)
      spls=f1
      go to 1000
c......
  400 if(l.ne.0) go to 1000
      spls=0.5d0*fgam(3-n)
      go to 1000
c......
  500 if(k.eq.4) go to 600
      if(cm.lt.1.d-8) go to 550
      call intga(gc,13)
      spls=f1
      go to 1000
  550 spls=0.5d0*fgam(3-n)*dexp(-pt)
      go to 1000
c......
  600 spls=0.5d0*fgam(3-n)
c......
 1000 continue
      return
      end
