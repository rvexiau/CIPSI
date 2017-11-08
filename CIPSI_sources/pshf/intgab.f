      subroutine intgab
      implicit real*8(a-h,o-z)
      common/uncp/ax,ay,az,bx,by,bz,cx,cy,cz,ga,gb,gc,am,bm,cm,
     1 pt,alfa,alfb,alfc,ctta,r,f1,f2,f3,f4,f5,f6,f7,f8,f9,xbs,bs(10),
     2 xg(64),hg(64),rxg(40),h1(40),h23(40),h456(40),h78(40),h9(40),
     3 fgam(13),rab,l,n,k,inum,nbs
      dimension f(9),index(6),tab(6)
      equivalence (f(1),f1)
      data tab/0.0,16.0,32.0,32.0,48.0,64.0/
      data np1,np2/40,64/
      data c/0.443113462726379d0/
      data index/1,2,4,3,6,9/
c...................
      if(n.eq.2) go to 300
      if(n.eq.-2) go to 200
      if(n.ne.0) go to 400
c......  n=0
  100 a=ga
      b=gb
      s=a+b
      p=a*b
      d=0.25d0*s*s-pt
      d=dexp(d)*c
      xbs=0.5d0*p
      go to (110,120,130,140,150,160),inum
  110 nbs=l
      call dfbsm1
      f1=d*bs(nbs+1)
      go to 1000
  120 nbs=l+1
      call dfbsm1
      f1=d*bs(nbs)
      f2=(a*bs(nbs)+b*bs(nbs+1))*d*0.5d0
      go to 1000
  130 nbs=l+2
      call dfbsm1
      d5=d*0.5d0
      f1=d*bs(nbs-1)
      f2=(a*bs(nbs-1)+b*bs(nbs))*d5
      f3=(b*bs(nbs-1)+a*bs(nbs))*d5
      f4=((bs(nbs-1)+bs(nbs+1))*p+(a*a+b*b+4.0d0)*bs(nbs))*d5*0.5d0
      go to 1000
  140 nbs=l+2
      call dfbsm1
      d5=d*0.5d0
      f1=d*bs(nbs-1)
      f2=(a*bs(nbs-1)+b*bs(nbs))*d5
      f3= (a*a*bs(nbs-1)+(p+p)*bs(nbs)+b*b*bs(nbs+1))*d5*0.5d0
      go to 1000
  150 nbs=l+3
      call dfbsm1
      d5=d*0.5d0
      d25=d5*0.5d0
      a2=a*a
      b2=b*b
      ab=a2+b2+4.0d0
      ab4=ab+4.0d0
      f1=d*bs(nbs-2)
      f2=(a*bs(nbs-2)+b*bs(nbs-1))*d5
      f3=(b*bs(nbs-2)+a*bs(nbs-1))*d5
      f4=(a2*bs(nbs-2)+(p+p)*bs(nbs-1)+b2*bs(nbs))*d25
      f5=((bs(nbs-2)+bs(nbs))*p+ab*bs(nbs-1))*d25
      f6=((p*bs(nbs-2)+(ab4+b2)*bs(nbs-1))*a
     1 +((ab4+a2)*bs(nbs)+p*bs(nbs+1))*b)*d25*0.5d0
      go to 1000
  160 nbs=l+4
      call dfbsm1
      d5=d*0.5d0
      d25=d5*0.5d0
      d125=d25*0.5d0
      a2=a*a
      b2=b*b
      ab2=a2+b2+2.0d0
      ab4=ab2+2.0d0
      ab8=ab4+2.0d0+2.0d0
      p2=p*p
      sp=p+p
      f1=d*bs(nbs-3)
      f2=(a*bs(nbs-3)+b*bs(nbs-2))*d5
      f3=(b*bs(nbs-3)+a*bs(nbs-2))*d5
      f4=(a2*bs(nbs-3)+sp*bs(nbs-2)+b2*bs(nbs-1))*d25
      f5=((bs(nbs-3)+bs(nbs-1))*p+ab4*bs(nbs-2))*d25
      f6=(b2*bs(nbs-3)+sp*bs(nbs-2)+a2*bs(nbs-1))*d25
      f7=((p*bs(nbs-3)+(ab8+b2)*bs(nbs-2))*a
     1 +((ab8+a2)*bs(nbs-1)+p*bs(nbs))*b)*d125
      f8=((p*bs(nbs-3)+(ab8+a2)*bs(nbs-2))*b
     1 +((ab8+b2)*bs(nbs-1)+p*bs(nbs))*a)*d125
      f9=((bs(nbs-3)+bs(nbs+1))*p2+(bs(nbs-2)+bs(nbs))*ab8*sp
     1 +((4.0d0*ab2+p2)*4.0d0+a2*a2+b2*b2)*bs(nbs-1))*d125*0.5d0
      go to 1000
c......  n=-2
  200 a=ga
      b=gb
      s=a+b
      p=a*b
      d=0.25d0*s*s-pt
      d=dexp(d)*c*0.25d0
      xl=dble(4*l)
      xbs=0.5d0*p
      go to (210,220,230,240,250,260),inum
  210 nbs=l+1
      call dfbsm1
      f1=((a*a+b*b+xl+6.0d0)*bs(nbs)+(p+p)*bs(nbs+1))*d
      go to 1000
  220 nbs=l+2
      call dfbsm1
      a2=a*a
      b2=b*b
      abl=a2+b2+xl
      sp=p+p
      f1=((abl+6.0d0)*bs(nbs-1)+sp*bs(nbs))*d
      f2=((abl+10.0d0)*a*bs(nbs-1)
     1 +((abl+a2+a2+14.0d0)*bs(nbs)+sp*bs(nbs+1))*b)*d*0.5d0
      go to 1000
  230 nbs=l+3
      call dfbsm1
      a2=a*a
      b2=b*b
      ab=a2+b2
      abl=ab+xl
      sp=p+p
      sp2=sp*p
      d5=d*0.5d0
      f1=((abl+6.0d0)*bs(nbs-2)+sp*bs(nbs-1))*d
      f2=((abl+10.0d0)*a*bs(nbs-2)+((abl+a2+a2+14.0d0)*bs(nbs-1)
     1 +sp*bs(nbs))*b)*d5
      f3=((abl+10.0d0)*b*bs(nbs-2)+((abl+b2+b2+14.0d0)*bs(nbs-1)
     1 +sp*bs(nbs))*a)*d5
      f4=(((abl+14.0d0)*bs(nbs-2)+(3.0d0*ab+xl+30.0d0)*bs(nbs))*p
     1 +((ab+4.0d0)*(abl+4.0d0)+14.0d0*ab+sp2+40.0d0)*bs(nbs-1)
     2 +sp2*bs(nbs+1))*d5*0.5d0
      go to 1000
  240 nbs=l+3
      call dfbsm1
      a2=a*a
      b2=b*b
      abl=a2+b2+xl
      sa2=a2+a2
      d5=d*0.5d0
      sp=p+p
      f1=((abl+6.0d0)*bs(nbs-2)+sp*bs(nbs-1))*d
      f2=((abl+10.0d0)*a*bs(nbs-2)+
     1  ((abl+sa2+14.0d0)*bs(nbs-1)+sp*bs(nbs))*b)*d5
      f3=((abl+14.0d0)*a2*bs(nbs-2)+(abl+a2+18.0d0)*sp*bs(nbs-1)
     1 +((abl+sa2+sa2+22.0d0)*bs(nbs)+sp*bs(nbs+1))*b2)*d5*0.5d0
      go to 1000
  250 nbs=l+4
      call dfbsm1
      a2=a*a
      b2=b*b
      ab=a2+b2
      abl=ab+xl
      sa2=a2+a2
      sb2=b2+b2
      sp=p+p
      sp2=p*sp
      d5=d*0.5d0
      d25=d5*0.5d0
      f1=((abl+6.0d0)*bs(nbs-3)+sp*bs(nbs-2))*d
      f2=((abl+10.0d0)*a*bs(nbs-3)
     1 +((abl+14.0d0+sa2)*bs(nbs-2)+sp*bs(nbs-1))*b)*d5
      f3=((abl+10.0d0)*b*bs(nbs-3)
     1 +((abl+14.0d0+sb2)*bs(nbs-2)+sp*bs(nbs-1))*a)*d5
      f4=((abl+14.0d0)*a2*bs(nbs-3)+(abl+a2+18.0d0)*sp*bs(nbs-2)
     1 +((abl+22.0d0+sa2+sa2)*bs(nbs-1)+sp*bs(nbs))*b2)*d25
      f5=(((abl+14.0d0)*bs(nbs-3)+(3.0d0*ab+xl+30.0d0)*bs(nbs-1))*p
     1 +((ab+4.0d0)*(abl+4.0d0)+14.0d0*ab+sp2+40.0d0)*bs(nbs-2)
     2 +sp2*bs(nbs))*d25
      f6=(((abl+18.0d0)*p*bs(nbs-3)+((ab+b2+8.0d0)*(abl+4.0d0)+sp2
     1 +(sb2+a2)*18.0d0+112.0d0)*bs(nbs-2))*a
     2+(((ab+a2+8.0d0)*(abl+8.0d0)+(ab+b2+26.0d0)*sa2+18.0d0*b2+112.0d0)
     3 *bs(nbs-1)+(3.0d0*ab+sa2+xl+46.0d0)*p*bs(nbs)+sp2*bs(nbs+1))*b)*
     4 d25*0.5d0
      go to 1000
  260 nbs=l+5
      call dfbsm1
      a2=a*a
      b2=b*b
      ab=a2+b2
      abl=ab+xl
      sa2=a2+a2
      sb2=b2+b2
      aabb=a2*a2+b2*b2
      sp=p+p
      sp2=sp*p
      ssp2=sp2+sp2
      d5=d*0.5d0
      d25=d5*0.5d0
      d125=d25*0.5d0
      f1=((abl+6.0d0)*bs(nbs-4)+sp*bs(nbs-3))*d
      f2=((abl+10.0d0)*a*bs(nbs-4)+((abl+14.0d0+sa2)*bs(nbs-3)
     1 +sp*bs(nbs-2))*b)*d5
      f3=((abl+10.0d0)*b*bs(nbs-4)+((abl+14.0d0+sb2)*bs(nbs-3)
     1 +sp*bs(nbs-2))*a)*d5
      f4=((abl+14.0d0)*a2*bs(nbs-4)+(abl+a2+18.0d0)*sp*bs(nbs-3)
     1 +((abl+22.0d0+sa2+sa2)*bs(nbs-2)+sp*bs(nbs-1))*b2)*d25
      f5=(((abl+14.0d0)*bs(nbs-4)+(3.0d0*ab+xl+30.0d0)*bs(nbs-2))*p
     1 +((ab+4.0d0)*(abl+4.0d0)+14.0d0*ab+sp2+40.0d0)*bs(nbs-3)
     2 +sp2*bs(nbs-1))*d25
      f6=((abl+14.0d0)*b2*bs(nbs-4)+(abl+18.0d0+b2)*sp*bs(nbs-3)
     1 +((abl+22.0d0+sb2+sb2)*bs(nbs-2)+sp*bs(nbs-1))*a2)*d25
      f7=(((abl+18.0d0)*p*bs(nbs-4)+((ab+b2+8.0d0)*(abl+4.0d0)+sp2
     1 +(sb2+a2)*18.0d0+112.0d0)*bs(nbs-3))*a+(((ab+a2+8.0d0)*
     3(abl+8.0d0)+(ab+b2+26.0d0)*sa2+18.0d0*b2+112.0d0)*bs(nbs-2)
     4 +(3.0d0*ab+sa2+xl+46.0d0)*p*bs(nbs-1)+sp2*bs(nbs))*b)*d125
      f8=(((abl+18.0d0)*p*bs(nbs-4)+((ab+a2+8.0d0)*(abl+4.0d0)+sp2
     1 +(sa2+b2)*18.0d0+112.0d0)*bs(nbs-3))*b+(((ab+b2+8.0d0)*
     3(abl+8.0d0)+(ab+a2+26.0d0)*sb2+18.0d0*a2+112.0d0)*bs(nbs-2)
     4 +(3.0d0*ab+sb2+xl+46.0d0)*p*bs(nbs-1)+sp2*bs(nbs))*a)*d125
      f9=((((abl+22.0d0)*bs(nbs-4)+(5.0d0*ab+xl+70.0d0)*bs(nbs))*p
     1 +sp2*bs(nbs+1))*p+(((ab+8.0d0)*(abl+4.0d0)+p*p+22.0d0*ab+144.0d0)
     2 *bs(nbs-3)+((ab+8.0d0)*(abl+12.0d0)+aabb+ssp2+38.0d0*ab+176.0d0)
     4 *bs(nbs-1))*sp+((abl+8.0d0)*(aabb+ssp2+16.0d0*ab+32.0d0)+ssp2*ab
     5 +60.0d0*sp2+22.0d0*aabb+288.0d0*ab+448.0d0)*bs(nbs-2))*d125*0.5d0
      go to 1000
c......  n=2
  300 imax=index(inum)
      do 305 i=1,imax
  305 f(i)=0.0d0
      do 340 ip=1,np1
      rx=rxg(ip)
      a=ga*rx
      b=gb*rx
      s=a+b
      p=a*b
      d=0.25d0*s*s-pt
      d=dexp(d)*c
      xbs=0.5d0*p
      go to (310,315,320,325,330,335),inum
  310 nbs=l
      call dfbsm1
      f1=d*bs(nbs+1)*h1(ip)+f1
      go to 340
  315 nbs=l+1
      call dfbsm1
      f1=d*bs(nbs) *h1(ip)+f1
      f2=(a*bs(nbs)+b*bs(nbs+1))*d*0.5d0*h23(ip)+f2
      go to 340
  320 nbs=l+2
      call dfbsm1
      d5=d*0.5d0
      f1=d*bs(nbs-1)*h1(ip)+f1
      f2=(a*bs(nbs-1)+b*bs(nbs))*d5*h23(ip)+f2
      f3=(b*bs(nbs-1)+a*bs(nbs))*d5*h23(ip)+f3
      f4=((bs(nbs-1)+bs(nbs+1))*p+(a*a+b*b+4.0d0)*bs(nbs))*d5*0.5d0
     1 *h456(ip)+f4
      go to 340
  325 nbs=l+2
      call dfbsm1
      d5=d*0.5d0
      f1=d*bs(nbs-1)*h1(ip) +f1
      f2=(a*bs(nbs-1)+b*bs(nbs))*d5*h23(ip)+f2
      f3= (a*a*bs(nbs-1)+(p+p)*bs(nbs)+b*b*bs(nbs+1))*d5*0.5d0*h456(ip)
     1  +f3
      go to 340
  330 nbs=l+3
      call dfbsm1
      d5=d*0.5d0
      d25=d5*0.5d0
      a2=a*a
      b2=b*b
      ab=a2+b2+4.0d0
      ab4=ab+4.0d0
      f1=d*bs(nbs-2)*h1(ip) +f1
      f2=(a*bs(nbs-2)+b*bs(nbs-1))*d5  *h23(ip) +f2
      f3=(b*bs(nbs-2)+a*bs(nbs-1))*d5  *h23(ip) +f3
      f4=(a2*bs(nbs-2)+(p+p)*bs(nbs-1)+b2*bs(nbs))*d25 *h456(ip) +f4
      f5=((bs(nbs-2)+bs(nbs))*p+ab*bs(nbs-1))*d25  *h456(ip)   +f5
      f6=((p*bs(nbs-2)+(ab4+b2)*bs(nbs-1))*a
     1 +((ab4+a2)*bs(nbs)+p*bs(nbs+1))*b)*d25*0.5d0  *h78(ip) +f6
      go to 340
  335 nbs=l+4
      call dfbsm1
      d5=d*0.5d0
      d25=d5*0.5d0
      d125=d25*0.5d0
      a2=a*a
      b2=b*b
      ab2=a2+b2+2.0d0
      ab4=ab2+2.0d0
      ab8=ab4+2.0d0+2.0d0
      p2=p*p
      sp=p+p
      f1=d*bs(nbs-3) *h1(ip) +f1
      f2=(a*bs(nbs-3)+b*bs(nbs-2))*d5  *h23(ip) +f2
      f3=(b*bs(nbs-3)+a*bs(nbs-2))*d5 *h23(ip)  +f3
      f4=(a2*bs(nbs-3)+sp*bs(nbs-2)+b2*bs(nbs-1))*d25   *h456(ip)+f4
      f5=((bs(nbs-3)+bs(nbs-1))*p+ab4*bs(nbs-2))*d25    *h456(ip) +f5
      f6=(b2*bs(nbs-3)+sp*bs(nbs-2)+a2*bs(nbs-1))*d25    *h456(ip)+f6
      f7=((p*bs(nbs-3)+(ab8+b2)*bs(nbs-2))*a
     1 +((ab8+a2)*bs(nbs-1)+p*bs(nbs))*b)*d125    *h78(ip) +f7
      f8=((p*bs(nbs-3)+(ab8+a2)*bs(nbs-2))*b
     1 +((ab8+b2)*bs(nbs-1)+p*bs(nbs))*a)*d125    *h78(ip) +f8
      f9=((bs(nbs-3)+bs(nbs+1))*p2+(bs(nbs-2)+bs(nbs))*ab8*sp
     1 +((4.0d0*ab2+p2)*4.0d0+a2*a2+b2*b2)*bs(nbs-1))*d125*0.5d0*h9(ip)
     2  +f9
  340 continue
      do 345 i=1,imax
  345 f(i)=f(i)*0.5d0
      go to 1000
c......  general case
  400 s=ga+gb
      dmin=8.0d0*dble(2*l+2-n)
      dmax=dmin+tab(inum)
      rmin=dsqrt(s*s+dmin)
      rmax=dsqrt(s*s+dmax)
      rmin=(rmin+s)*0.25d0
      rmax=(rmax+s)*0.25d0
      if(rmin.gt.6.0d0) go to 405
      bpa=(rmax+6.0d0)*0.5d0
      bma=bpa
      go to 410
  405 bpa=(rmax+rmin)*0.5d0
      bma=(rmax-rmin+12.0d0)*0.5d0
  410 imax=index(inum)
      do 415 i=1,imax
  415 f(i)=0.0d0
      do 450 i=1,np2
      x=bma*xg(i)+bpa
      t=dexp((s-x)*x-pt)
      t=t*x**(2-n)
      t=t*hg(i)
      go to (420,425,430,435,440,445),inum
  420 nbs=l
      xbs=ga*x
      call dfbsm1
      a1=bs(nbs+1)
      xbs=gb*x
      call dfbsm1
      b1=bs(nbs+1)
      f1=f1+t*a1*b1
      go to 450
  425 nbs=l+1
      xbs=ga*x
      call dfbsm1
      a1=bs(nbs)
      a2=bs(nbs+1)
      nbs=l
      xbs=gb*x
      call dfbsm1
      b1=bs(nbs+1)
      f1=f1+t*a1*b1
      f2=f2+t*a2*b1*x
      go to 450
  430 nbs=l+1
      xbs=ga*x
      call dfbsm1
      a1=bs(nbs)
      a2=bs(nbs+1)
      xbs=gb*x
      call dfbsm1
      b1=bs(nbs)
      b2=bs(nbs+1)
      f1=f1+t*a1*b1
      t=t*x
      f2=f2+t*a2*b1
      f3=f3+t*a1*b2
      t=t*x
      f4=f4+t*a2*b2
      go to 450
  435 nbs=l+2
      xbs=ga*x
      call dfbsm1
      a1=bs(nbs-1)
      a2=bs(nbs)
      a3=bs(nbs+1)
      nbs=l
      xbs=gb*x
      call dfbsm1
      b1=bs(nbs+1)
      f1=f1+t*a1*b1
      f2=f2+t*a2*b1*x
      f3=f3+t*a3*b1*x*x
      go to 450
  440 nbs=l+2
      xbs=ga*x
      call dfbsm1
      a1=bs(nbs-1)
      a2=bs(nbs)
      a3=bs(nbs+1)
      nbs=l+1
      xbs=gb*x
      call dfbsm1
      b1=bs(nbs)
      b2=bs(nbs+1)
      f1=f1+t*a1*b1
      t=t*x
      f2=f2+t*a2*b1
      f3=f3+t*a1*b2
      t=t*x
      f4=f4+t*a3*b1
      f5=f5+t*a2*b2
      t=t*x
      f6=f6+t*a3*b2
      go to 450
  445 nbs=l+2
      xbs=ga*x
      call dfbsm1
      a1=bs(nbs-1)
      a2=bs(nbs)
      a3=bs(nbs+1)
      xbs=gb*x
      call dfbsm1
      b1=bs(nbs-1)
      b2=bs(nbs)
      b3=bs(nbs+1)
      f1=f1+t*a1*b1
      t=t*x
      f2=f2+t*a2*b1
      f3=f3+t*a1*b2
      t=t*x
      f4=f4+t*a3*b1
      f5=f5+t*a2*b2
      f6=f6+t*a1*b3
      t=t*x
      f7=f7+t*a3*b2
      f8=f8+t*a2*b3
      t=t*x
      f9=f9+t*a3*b3
  450 continue
      do 455 i=1,imax
  455 f(i)=f(i)*bma
c...................
 1000 return
      end
