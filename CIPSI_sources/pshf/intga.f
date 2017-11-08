      subroutine intga(a,icas)
      implicit real*8(a-h,o-z)
      common/uncp/ax,ay,az,bx,by,bz,cx,cy,cz,ga,gb,gc,am,bm,cm,
     1 pt,alfa,alfb,alfc,ctta,r,f1,f2,f3,f4,f5,f6,f7,f8,f9,xbs,bs(10),
     2 xg(64),hg(64),rxg(40),h1(40),h23(40),h456(40),h78(40),h9(40),
     3 fgam(13),rab,l,n,k,inum,nbs
      dimension f(9),index(23),dml(23),dmu(23)
      equivalence (f1,f(1))
      data index/1,1,2,2,1,1,3,2,2,3,3,3,1,2,1,3,2,3,2,4,3,5,4/
      data np1,np2/40,64/
      data dml/0.0,16.0,0.0,16.0,16.0,32.0,0.0,16.0,32.0,16.0,16.0,
     1 32.,0.,0.,16.,0.,16.,0.,16.,0.,16.,0.,16./
      data dmu/0.0,0.0,16.0,16.0,0.0,0.0,32.0,16.0,16.0,32.0,32.0,32.0,
     1 0.,16.,16.,32.,32.,32.,32.,48.,48.,64.,64./
      data c0/0.4431134627263797d0/
      data c1/0.2215567313631895d0/
      data c2/0.1107783656815948d0/
      data c3/0.5538918284079738d-1/
      data c4/0.2769459142039869d-1/
      data c5/0.1384729571019934d-1/
      data c6/0.6923647855099672d-2/
c......
      if(n.eq.2) go to 300
      if(n.eq.-2) go to 200
      if(n.ne.0) go to 400
c......  n.eq.0
  100 a2=a*a
      d=dexp(0.25d0*a2-pt)
      go to (102,104,106,108,110,112,114,116,118,120,122,124,126,128,
     1 130,132,134,136,138,140,142,144,146),icas
  102 f1=c0*d
      go to 1000
  104 f1=c1*d*a
      go to 1000
  106 f1=c0*d
      f2=c1*d*a
      go to 1000
  108 f1=c1*d*a
      f2=c2*d*a2
      go to 1000
  110 f1=(a2+6.0d0)*c2*d
      go to 1000
  112 f1=c2*d*a2
      go to 1000
  114 f1=c0*d
      f2=c1*d*a
      f3=c2*d*a2
      go to 1000
  116 f1=(a2+6.0d0)*c2*d
      f2=(a2+10.0d0)*c3*d*a
      go to 1000
  118 f1=c2*d*a2
      f2=c3*d*a2*a
      go to 1000
  120 f1=c1*d*a
      f2=c2*d*a2
      f3=c3*d*a2*a
      go to 1000
  122 f1=(a2+6.0d0)*c2*d
      f2=(a2+10.0d0)*c3*d*a
      f3=(a2+14.0d0)*c4*d*a2
      go to 1000
  124 f1=c2*d*a2
      f2=c3*d*a2*a
      f3=c4*d*a2*a2
      go to 1000
  126 f1=c0*d
      go to 1000
  128 f1=c0*d
      f2=c1*d*a
      go to 1000
  130 f1=c1*d*a
      go to 1000
  132 f1=c0*d
      f2=c1*d*a
      f3=c2*d*a2
      go to 1000
  134 f1=c1*d*a
      f2=c2*d*a2
      go to 1000
  136 f1=c0*d
      f2=c1*d*a
      f3=c2*d*a2
      go to 1000
  138 f1=c1*d*a
      f2=c2*d*a2
      go to 1000
  140 f1=c0*d
      f2=c1*d*a
      f3=c2*d*a2
      f4=c3*d*a2*a
      go to 1000
  142 f1=c1*d*a
      f2=c2*d*a2
      f3=c3*d*a2*a
      go to 1000
  144 f1=c0*d
      f2=c1*d*a
      f3=c2*d*a2
      f4=c3*d*a2*a
      f5=c4*d*a2*a2
      go to 1000
  146 f1=c1*d*a
      f2=c2*d*a2
      f3=c3*d*a2*a
      f4=c4*d*a2*a2
      go to 1000
c......  n.eq. -2
  200 a2=a*a
      d=dexp(0.25d0*a2-pt)
      go to (202,204,206,208,210,212,214,216,218,220,222,224,226,
     1 228,230,232,234,236,238,240,242,244,246),icas
  202 f1=(a2+6.0d0)*c2*d
      go to 1000
  204 f1=(a2+10.0d0)*c3*d*a
      go to 1000
  206 f1=(a2+6.0d0)*c2*d
      f2=(a2+10.0d0)*c3*d*a
      go to 1000
  208 f1=(a2+10.0d0)*c3*d*a
      f2=(a2+14.0d0)*c4*d*a2
      go to 1000
  210 f1=((a2+20.0d0)*a2+60.0d0)*c4*d
      go to 1000
  212 f1=(a2+14.0d0)*c4*d*a2
      go to 1000
  214 f1=(a2+6.0d0)*c2*d
      f2=(a2+10.0d0)*c3*d*a
      f3=(a2+14.0d0)*c4*d*a2
      go to 1000
  216 f1=((a2+20.0d0)*a2+60.0d0)*c4*d
      f2=((a2+28.0d0)*a2+140.0d0)*c5*d*a
      go to 1000
  218 f1=(a2+14.0d0)*c4*d*a2
      f2=(a2+18.0d0)*c5*d*a2*a
      go to 1000
  220 f1=(a2+10.0d0)*c3*d*a
      f2=(a2+14.0d0)*c4*d*a2
      f3=(a2+18.0d0)*c5*d*a2*a
      go to 1000
  222 f1=((a2+20.0d0)*a2+60.0d0)*c4*d
      f2=((a2+28.0d0)*a2+140.0d0)*c5*d*a
      f3=((a2+36.0d0)*a2+252.0d0)*c6*d*a2
      go to 1000
  224 f1=(a2+14.0d0)*c4*d*a2
      f2=(a2+18.0d0)*c5*d*a2*a
      f3=(a2+22.0d0)*c6*d*a2*a2
      go to 1000
  226 f1=(a2+6.0d0)*c2*d
      go to 1000
  228 f1=(a2+6.0d0)*c2*d
      f2=(a2+10.0d0)*c3*d*a
      go to 1000
  230 f1=(a2+10.0d0)*c3*d*a
      go to 1000
  232 f1=(a2+6.0d0)*c2*d
      f2=(a2+10.0d0)*c3*d*a
      f3=(a2+14.0d0)*c4*d*a2
      go to 1000
  234 f1=(a2+10.0d0)*c3*d*a
      f2=(a2+14.0d0)*c4*d*a2
      go to 1000
  236 f1=(a2+6.0d0)*c2*d
      f2=(a2+10.0d0)*c3*d*a
      f3=(a2+14.0d0)*c4*d*a2
      go to 1000
  238 f1=(a2+10.0d0)*c3*d*a
      f2=(a2+14.0d0)*c4*d*a2
      go to 1000
  240 f1=(a2+6.0d0)*c2*d
      f2=(a2+10.0d0)*c3*d*a
      f3=(a2+14.0d0)*c4*d*a2
      f4=(a2+18.0d0)*c5*d*a2*a
      go to 1000
  242 f1=(a2+10.0d0)*c3*d*a
      f2=(a2+14.0d0)*c4*d*a2
      f3=(a2+18.0d0)*c5*d*a2*a
      go to 1000
  244 f1=(a2+6.0d0)*c2*d
      f2=(a2+10.0d0)*c3*d*a
      f3=(a2+14.0d0)*c4*d*a2
      f4=(a2+18.0d0)*c5*d*a2*a
      f5=(a2+22.0d0)*c6*d*a2*a2
      go to 1000
  246 f1=(a2+10.0d0)*c3*d*a
      f2=(a2+14.0d0)*c4*d*a2
      f3=(a2+18.0d0)*c5*d*a2*a
      f4=(a2+22.0d0)*c6*d*a2*a2
      go to 1000
c......  n.eq.2
  300 go to (308,308,308,308,302,308,308,304,308,308,306,308,
     1 308,308,308,308,308,308,308,308,308,308,308),icas
  302 f1=c0*dexp(0.25d0*a*a-pt)
      go to 1000
  304 d=dexp(0.25d0*a*a-pt)
      f1=c0*d
      f2=c1*d*a
      go to 1000
  306 a2=a*a
      d=dexp(0.25d0*a2-pt)
      f1=c0*d
      f2=c1*d*a
      f3=c2*d*a2
      go to 1000
  308 imax=index(icas)
      do 310 i=1,imax
  310 f(i)=0.0d0
      do 350 i=1,np1
      x=a*rxg(i)
      x2=x*x
      d=dexp(0.25d0*x2-pt)
      go to (312,314,316,318,1000,320,322,1000,324,326,1000,328,
     1 330,332,334,336,338,340,342,344,346,348,349),icas
  312 f1=f1+c0*d*h1(i)
      go to 350
  314 f1=f1+c1*d*x*h23(i)
      go to 350
  316 f1=f1+c0*d*h1(i)
      f2=f2+c1*d*x*h23(i)
      go to 350
  318 f1=f1+c1*d*x*h23(i)
      f2=f2+c2*d*x2*h456(i)
      go to 350
  320 f1=f1+c2*d*x2*h456(i)
      go to 350
  322 f1=f1+c0*d*h1(i)
      f2=f2+c1*d*x*h23(i)
      f3=f3+c2*d*x2*h456(i)
      go to 350
  324 f1=f1+c2*d*x2*h456(i)
      f2=f2+c3*d*x2*x*h78(i)
      go to 350
  326 f1=f1+c1*d*x*h23(i)
      f2=f2+c2*d*x2*h456(i)
      f3=f3+c3*d*x2*x*h78(i)
      go to 350
  328 f1=f1+c2*d*x2*h456(i)
      f2=f2+c3*d*x2*x*h78(i)
      f3=f3+c4*d*x2*x2*h9(i)
      go to 350
  330 f1=f1+c0*d*h1(i)
      go to 350
  332 f1=f1+c0*d*h1(i)
      f2=f2+c1*d*x*h23(i)
      go to 350
  334 f1=f1+c1*d*x*h23(i)
      go to 350
  336 f1=f1+c0*d*h1(i)
      f2=f2+c1*d*x*h23(i)
      f3=f3+c2*d*x2*h456(i)
      go to 350
  338 f1=f1+c1*d*x*h23(i)
      f2=f2+c2*d*x2*h456(i)
      go to 350
  340 f1=f1+c0*d*h1(i)
      f2=f2+c1*d*x*h23(i)
      f3=f3+c2*d*x2*h456(i)
      go to 350
  342 f1=f1+c1*d*x*h23(i)
      f2=f2+c2*d*x2*h456(i)
      go to 350
  344 f1=f1+c0*d*h1(i)
      f2=f2+c1*d*x*h23(i)
      f3=f3+c2*d*x2*h456(i)
      f4=f4+c3*d*x2*x*h78(i)
      go to 350
  346 f1=f1+c1*d*x*h23(i)
      f2=f2+c2*d*x2*h456(i)
      f3=f3+c3*d*x2*x*h78(i)
      go to 350
  348 f1=f1+c0*d*h1(i)
      f2=f2+c1*d*x*h23(i)
      f3=f3+c2*d*x2*h456(i)
      f4=f4+c3*d*x2*x*h78(i)
      f5=f5+c4*d*x2*x2*h9(i)
  349 f1=f1+c1*d*x*h23(i)
      f2=f2+c2*d*x2*h456(i)
      f3=f3+c3*d*x2*x*h78(i)
      f4=f4+c4*d*x2*x2*h9(i)
      go to 350
  350 continue
      do 360 i=1,imax
  360 f(i)=f(i)*0.5d0
      go to 1000
c......  general case
  400 dm=(2-n)*8
      dmin=dm+dml(icas)
      dmax=dmin+dmu(icas)
      rmin=dsqrt(a*a+dmin)
      rmax=dsqrt(a*a+dmax)
      rmin=(a+rmin)*0.25d0
      rmax=(a+rmax)*0.25d0
      if(rmin.gt.6.0d0) go to 405
      bpa=(rmax+6.0d0)*0.5d0
      bma=bpa
      go to 410
  405 bpa=(rmax+rmin)*0.5d0
      bma=(rmax-rmin+12.0d0)*0.5d0
  410 imax=index(icas)
      do 415 i=1,imax
  415 f(i)=0.0d0
      do 470 i=1,np2
      x=bma*xg(i)+bpa
      x2=x*x
      xbs=a*x
      t=dexp((a-x)*x-pt)
      t=t*x**(2-n)
      t=t*hg(i)
      go to (422,424,426,428,430,432,434,436,438,440,442,444,
     1 446,448,450,452,454,456,458,460,462,464,466),icas
  422 nbs=0
      call dfbsm1
      f1=f1+t*bs(nbs+1)
      go to 470
  424 nbs=1
      call dfbsm1
      f1=f1+t*bs(nbs+1)*x
      go to 470
  426 nbs=1
      call dfbsm1
      f1=f1+t*bs(nbs)
      f2=f2+t*bs(nbs+1)*x
      go to 470
  428 nbs=2
      call dfbsm1
      f1=f1+t*bs(nbs)*x
      f2=f2+t*bs(nbs+1)*x2
      go to 470
  430 nbs=0
      call dfbsm1
      f1=f1+t*bs(nbs+1)*x2
      go to 470
  432 nbs=2
      call dfbsm1
      f1=f1+t*bs(nbs+1)*x2
      go to 470
  434 nbs=2
      call dfbsm1
      f1=f1+t*bs(nbs-1)
      f2=f2+t*bs(nbs)*x
      f3=f3+t*bs(nbs+1)*x2
      go to 470
  436 nbs=1
      call dfbsm1
      f1=f1+t*bs(nbs)*x2
      f2=f2+t*bs(nbs+1)*x2*x
      go to 470
  438 nbs=3
      call dfbsm1
      f1=f1+t*bs(nbs)*x2
      f2=f2+t*bs(nbs+1)*x2*x
      go to 470
  440 nbs=3
      call dfbsm1
      f1=f1+t*bs(nbs-1)*x
      f2=f2+t*bs(nbs)*x2
      f3=f3+t*bs(nbs+1)*x2*x
      go to 470
  442 nbs=2
      call dfbsm1
      f1=f1+t*bs(nbs-1)*x2
      f2=f2+t*bs(nbs)*x2*x
      f3=f3+t*bs(nbs+1)*x2*x2
      go to 470
  444 nbs=4
      call dfbsm1
      f1=f1+t*bs(nbs-1)*x2
      f2=f2+t*bs(nbs)*x2*x
      f3=f3+t*bs(nbs+1)*x2*x2
      go to 470
  446 nbs=0
      call dfbsm1
      f1=f1+t*bs(nbs+1)
      go to 470
  448 nbs=1
      call dfbsm1
      f1=f1+t*bs(nbs)
      f2=f2+t*bs(nbs+1)*x
      go to 470
  450 nbs=1
      call dfbsm1
      f1=f1+t*bs(nbs+1)*x
      go to 470
  452 nbs=2
      call dfbsm1
      f1=f1+t*bs(nbs-1)
      f2=f2+t*bs(nbs)*x
      f3=f3+t*bs(nbs+1)*x2
      go to 470
  454 nbs=2
      call dfbsm1
      f1=f1+t*bs(nbs)*x
      f2=f2+t*bs(nbs+1)*x2
      go to 470
  456 nbs=2
      call dfbsm1
      f1=f1+t*bs(nbs-1)
      f2=f2+t*bs(nbs)*x
      f3=f3+t*bs(nbs+1)*x2
      go to 470
  458 nbs=2
      call dfbsm1
      f1=f1+t*bs(nbs)*x
      f2=f2+t*bs(nbs+1)*x2
      go to 470
  460 nbs=3
      call dfbsm1
      f1=f1+t*bs(nbs-2)
      f2=f2+t*bs(nbs-1)*x
      f3=f3+t*bs(nbs)*x2
      f4=f4+t*bs(nbs+1)*x2*x
      go to 470
  462 nbs=3
      call dfbsm1
      f1=f1+t*bs(nbs-1)*x
      f2=f2+t*bs(nbs)*x2
      f3=f3+t*bs(nbs+1)*x2*x
      go to 470
  464 nbs=4
      call dfbsm1
      f1=f1+t*bs(nbs-3)
      f2=f2+t*bs(nbs-2)*x
      f3=f3+t*bs(nbs-1)*x2
      f4=f4+t*bs(nbs)*x2*x
      f5=f5+t*bs(nbs+1)*x2*x2
      go to 470
  466 nbs=4
      call dfbsm1
      f1=f1+t*bs(nbs-2)*x
      f2=f2+t*bs(nbs-1)*x2
      f3=f3+t*bs(nbs)*x2*x
      f4=f4+t*bs(nbs+1)*x2*x2
  470 continue
      do 480 i=1,imax
  480 f(i)=f(i)*bma
c......
 1000 return
      end
