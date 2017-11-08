      subroutine pgdplp
       implicit real*8 (a-h,o-z)
      common/uncp/a(3),b(3),c(3),ga,gb,gc,am,bm,cm,
     1 pt,alfa,alfb,alfc,ctta,r,f1,f2,f3,f4,f5,f6,f7,f8,f9,xbs,bs(10),
     2 xg(64),hg(64),rxg(40),h1(40),h23(40),h456(40),h78(40),h9(40),
     3 fgam(13),rab,l,n,k,inum,nbs
      common/pmat/psdint(225)
      dimension id(6),jd(6)
      data id/1,2,3,1,1,2/
      data jd/1,2,3,2,3,3/
      do 5 i=1,18
    5 psdint(i)=0.d0
c......
c......
      if(l.lt.0) go to 500
      go to (100,200,300,1000),k
  100 call intgab
      ij=0
      do 150 m=1,6
      i=id(m)
      j=jd(m)
      do 150 jb=1,3
      ij=ij+1
      aab=a(i)*a(j)*b(jb)
      aaa=a(i)*a(j)*a(jb)
      abm=am*bm
      am2=am*am
      bm2=bm*bm
      al=2.d0*alfa*am
      alm=2.d0*alfa*abm
      if(l-1) 110,120,130
  110 t10=-aab
      t11=2.d0*aab/am
      if(i.eq.j) t11=t11-b(jb)/al
      t12=aab/bm
      t13=-aab/am2
      t14=-2.d0*aab/abm
      if(i.eq.j) t14=t14+b(jb)/alm
      t15=aab/(am*abm)
      go to 140
  120 ab=a(i)*b(j)+a(j)*b(i)
      bet=2.d0*alfb*ga
      w=aab*ctta
      dam=am*gb
      dbm=bm*ga
      d=b(jb)*ab/dbm
      e=aaa/dam
      t10=-w+e+d
      if(i.eq.jb) t10=t10-a(j)/(bet*bm)
      if(j.eq.jb) t10=t10-a(i)/(bet*bm)
      t11=(2.d0*(w-e)-d)/am
      if(i.eq.j) t11=t11+(a(jb)/dam-b(jb)*ctta)/al
      if(i.eq.jb) t11=t11+a(j)/(bet*abm)
      if(j.eq.jb) t11=t11+a(i)/(bet*abm)
      t12=(w-d)/bm
      t13=(e-w)/am2
      t14=(d-2.d0*w)/abm
      if(i.eq.j) t14=t14+b(jb)*ctta/alm
      t15=w/(abm*am)
      go to 140
  130 pl=1.5d0*ctta*ctta-0.5d0
      ab=a(i)*b(j)+a(j)*b(i)
      bbb=b(i)*b(j)*b(jb)
      ua=1.d0/(am*ga)
      ub=1.d0/(bm*gb)
      ca=3.d0*ctta*a(jb)
      cb=3.d0*ctta*b(jb)
      uam=1.d0/(am*gb)
      ubm=1.d0/(bm*ga)
      t10=-aab*(pl+2.d0*ua+ub-4.d0*ua*ub)
      t10=t10+3.d0*aaa*ctta*uam-3.d0*bbb*ubm*ubm
      t10=t10+cb*ab*ubm          -3.d0*a(jb)*ab*ua*ub
      if(i.eq.j) t10=t10+b(jb)*(1.d0-2.d0*ub)/(al*ga)
      if(i.eq.jb) t10=t10+3.d0*(b(j)*ubm-a(j)*ctta)/(r*ga*gb)
      if(j.eq.jb) t10=t10+3.d0*(b(i)*ubm-a(i)*ctta)/(r*ga*gb)
      t11=(2.d0*aab*(pl+ua+ub-2.d0*ua*ub)-6.d0*aaa*ctta*uam)/am
      t11=t11+ab*ua*(3.d0*a(jb)*uam-cb)/bm
      if(i.eq.j) t11=t11+(ca*uam-b(jb)*(pl+ub))/al
      if(i.eq.jb) t11=t11+3.d0*a(j)*ctta*ua/(r*gb)
      if(j.eq.jb) t11=t11+3.d0*a(i)*ctta*ua/(r*gb)
      t12=aab*(pl+2.d0*ua)/bm+(3.d0*bbb*ubm-ab*cb)*ubm/bm
      if(i.eq.j) t12=t12-b(jb)*ubm/al
      t13=(3.d0*aaa*ctta*uam-aab*(pl+ub))/am2
      t14=(ab*cb*ubm-2.d0*aab*(pl+ua))/abm
      if(i.eq.j) t14=t14+b(jb)*pl/(al*bm)
      t15=aab*pl/(am*abm)
  140 continue
      dplp=t10*f1+t11*f2+t12*f3+t13*f4+t14*f5+t15*f6
      psdint(ij)=dplp
  150 continue
      go to 1000
c......
  200 if(l-1) 210,1000,220
  210 call intga(gb,8)
      fy=(-f1+f2/bm)/(3.d0*r)
      ij=0
      do 215 m=1,3
      do 215 jb=1,3
      ij=ij+1
      psdint(ij)=b(jb)*fy
  215 continue
      go to 1000
  220 bm2=bm*bm
      call intga(gb,9)
      ij=0
      do 225 m=1,6
      i=id(m)
      j=jd(m)
      do 225 jb=1,3
      ij=ij+1
      bbbm=b(i)*b(j)*b(jb)/(5.d0*bm2)
      t1=-bbbm
      if(i.eq.j) t1=t1+b(jb)*(1.d0-2.d0/(gb*bm))/(15.d0*r)
      if(i.eq.jb) t1=t1+b(j)/(5.d0*r*gb*bm)
      if(j.eq.jb) t1=t1+b(i)/(5.d0*r*gb*bm)
      t2=bbbm/bm
      if(i.eq.j) t2=t2-b(jb)/(15.d0*r*bm)
      dplp=t1*f1+t2*f2
      psdint(ij)=dplp
  225 continue
      go to 1000
c......
  300 if(l.ne.1) go to 1000
      call intga(ga,10)
      ij=0
      do 310 m=1,6
      i=id(m)
      j=jd(m)
      do 310 jb=1,3
      ij=ij+1
      aaa=a(i)*a(j)*a(jb)
      am2=am*am
      al=2.0d0*alfa*am
      t1=aaa/am
      if(i.eq.jb) t1=t1-a(j)/(ga*r)
      if(j.eq.jb) t1=t1-a(i)/(ga*r)
      t2=-2.d0*aaa/am2
      if(i.eq.j) t2=t2+a(jb)/(al*am)
      if(i.eq.jb) t2=t2+a(j)/(ga*am*r)
      if(j.eq.jb) t2=t2+a(i)/(ga*am*r)
      t3=aaa/(am2*am)
      dplp=(t1*f1+t2*f2+t3*f3)*0.33333333333333333d0
      psdint(ij)=dplp
  310 continue
      go to 1000
c......
  500 go to (600,700,800,1000),k
  600 if(cm.lt.1.d-8) go to 650
      w=0.5d0/rab
      s=1.d0/cm
      call intga(gc,20)
      ij=0
      do 605 m=1,6
      i=id(m)
      j=jd(m)
      do 605 jb=1,3
      ij=ij+1
      cc=c(i)*c(j)
      caac=c(i)*a(j)+a(i)*c(j)
      t1=-a(i)*a(j)*b(jb)
      t2=a(i)*a(j)*c(jb)+caac*b(jb)
      t3=-cc*b(jb)-caac*c(jb)
      t4=cc*c(jb)
      if(i.ne.j) go to 610
      t2=t2-w*b(jb)
      t3=t3+w*c(jb)
  610 if(i.ne.jb) go to 620
      t2=t2-w*a(j)
      t3=t3+w*c(j)
  620 if(j.ne.jb) go to 630
      t2=t2-w*a(i)
      t3=t3+w*c(i)
  630 continue
      dplp=((t4*f4*s+t3*f3)*s+t2*f2)*s+t1*f1
      psdint(ij)=dplp
  605 continue
       go to 1000
c..
  650 ij=0
      do 655 m=1,6
      i=id(m)
      j=jd(m)
      do 655 jb=1,3
      ij=ij+1
      s=0.0d0
      if(i.eq.j) s=s+b(jb)
      if(i.eq.jb) s=s+a(j)
      if(j.eq.jb) s=s+a(i)
      dplp=-fgam(5-n)*s/(6.0d0*r)
      dplp=dplp-0.5d0*a(i)*a(j)*b(jb)*fgam(3-n)
      dplp=dplp*dexp(-pt)
      psdint(ij)=dplp
  655 continue
      go to 1000
c......
  700 if(cm.lt.1.d-8) go to 750
      w=0.5d0/rab
      s=1.d0/cm
      call intga(gc,21)
      ij=0
      do 705 m=1,6
      i=id(m)
      j=jd(m)
      do 705 jb=1,3
      ij=ij+1
      cc=c(i)*c(j)
      t1=0.0d0
      t2=-cc*b(jb)
      t3=cc*c(jb)
      if(i.ne.j) go to 710
      t1=t1-w*b(jb)
      t2=t2+w*c(jb)
  710 if(i.eq.jb) t2=t2+w*c(j)
      if(j.eq.jb) t2=t2+w*c(i)
      dplp=((t3*f3*s+t2*f2)*s+t1*f1)*s
      psdint(ij)=dplp
  705 continue
      go to 1000
  750 dplp=-fgam(5-n)/(6.0d0*r)
      dplp=dplp*dexp(-pt)
      ij=0
      do 755 m=1,3
      i=id(m)
      j=jd(m)
      do 755 jb=1,3
      ij=ij+1
      psdint(ij)=b(jb)*dplp
  755 continue
      go to 1000
c......
  800 if(cm.lt.1.d-8) go to 850
      w=0.5d0/rab
      s=1.d0/cm
      call intga(gc,21)
      ij=0
      do 805 m=1,6
      i=id(m)
      j=jd(m)
      do 805 jb=1,3
      ij=ij+1
      cc=c(i)*c(j)
      caac=c(i)*a(j)+a(i)*c(j)
      t1=a(i)*a(j)*c(jb)
      t2=-caac*c(jb)
      t3=cc*c(jb)
      if(i.eq.j) t2=t2+w*c(jb)
      if(i.ne.jb) go to 810
      t1=t1-w*a(j)
      t2=t2+w*c(j)
  810 if(j.ne.jb) go to 820
      t1=t1-w*a(i)
      t2=t2+w*c(i)
  820 continue
      dplp=((t3*f3*s+t2*f2)*s+t1*f1)*s
      psdint(ij)=dplp
  805 continue
      go to 1000
  850 dplp=-fgam(5-n)*dexp(-pt)/(6.d0*r)
      ij=0
      do 855 m=1,6
      i=id(m)
      j=jd(m)
      do 855 jb=1,3
      ij=ij+1
      s=0.d0
      if(i.eq.jb) s=s+a(j)
      if(j.eq.jb) s=s+a(i)
      psdint(ij)=s*dplp
  855 continue
c.....
 1000 continue
      return
      end
