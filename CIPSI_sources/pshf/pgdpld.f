      subroutine pgdpld
       implicit real*8 (a-h,o-z)
      logical*1 iandj
      common/uncp/a(3),b(3),c(3),ga,gb,gc,am,bm,cm,
     1 pt,alfa,alfb,alfc,ctta,r,f1,f2,f3,f4,f5,f6,f7,f8,f9,xbs,bs(10),
     2 xg(64),hg(64),rxg(40),h1(40),h23(40),h456(40),h78(40),h9(40),
     3 fgam(13),rab,l,n,k,inum,nbs
      common/pmat/ psdint(225),pspt(225),iandj
      dimension id(6),jd(6)
      data id/1,2,3,1,1,2/,jd/1,2,3,2,3,3/
      do 5 i=1,36
    5 psdint(i)=0.d0
      maxk=6
c......
      if(l.lt.0) go to 500
      go to (100,200,300,400),k
c......
  100 call intgab
      ij=0
      do 150 m=1,6
      i=id(m)
      j=jd(m)
      if(iandj)maxk=m
      do 150 mk=1,maxk
      ik=id(mk)
      im=jd(mk)
      ij=ij+1
      aa=a(i)*a(j)
      bb=b(ik)*b(im)
      aabb=aa*bb
      ala=2.d0*alfa*am
      alb=2.d0*alfb*bm
      abm=am*bm
      am2=am*am
      bm2=bm*bm
      if(l-1)110,120,130
  110 t1=aabb
      t2=-2.d0*aabb/am
      if(i.eq.j) t2=t2+bb/ala
      t3=-2.d0*aabb/bm
      if(ik.eq.im) t3=t3+aa/alb
      t4=aabb/am2
      t5=4.d0*aabb/abm
      if(i.eq.j) t5=t5-2.d0*bb/(ala*bm)
      if(i.eq.j.and.ik.eq.im) t5=t5+1.d0/(ala*alb)
      if(ik.eq.im) t5=t5-2.d0*aa/(alb*am)
      t6=aabb/bm2
      t7=-2.d0*aabb/(abm*am)
      if(ik.eq.im) t7=t7+aa/(am2*alb)
      t8=-2.d0*aabb/(abm*bm)
      if(i.eq.j) t8=t8+bb/(bm2*ala)
      t9=aabb/(abm*abm)
      go to 140
  120 w=aabb*ctta
      aibj=a(i)*b(j)
      ajbi=a(j)*b(i)
      akbm=a(ik)*b(im)
      ambk=a(im)*b(ik)
      aibm=a(i)*b(im)
      aibk=a(i)*b(ik)
      ajbm=a(j)*b(im)
      ajbk=a(j)*b(ik)
      rkm=akbm+ambk
      rij=aibj+ajbi
      raa=aa*rkm
      rbb=bb*rij
      dam=am*gb
      dbm=bm*ga
      gab=ga*gb
      alad=ala*dam
      albd=alb*dbm
      albe=albd*am
      alae=alad*bm
      t1=w-rbb/dbm-raa/dam
      if(i.eq.ik) t1=t1+ajbm/(gab*r)
      if(i.eq.im) t1=t1+ajbk/(gab*r)
      if(j.eq.ik) t1=t1+aibm/(gab*r)
      if(j.eq.im) t1=t1+aibk/(gab*r)
      t2=-2.d0*w/am+2.d0*raa/(dam*am)+rbb/(ga*abm)
      if(i.eq.j) t2=t2+(bb*ctta-rkm/dam)/ala
      if(i.eq.ik) t2=t2-ajbm/alad
      if(j.eq.ik) t2=t2-aibm/alad
      if(i.eq.im) t2=t2-ajbk/alad
      if(j.eq.im) t2=t2-aibk/alad
      t3=-2.d0*w/bm+2.d0*rbb/(dbm*bm)+raa/(gb*abm)
      if(ik.eq.im) t3=t3+(aa*ctta-rij/dbm)/alb
      if(i.eq.ik) t3=t3-ajbm/albd
      if(j.eq.ik) t3=t3-aibm/albd
      if(i.eq.im) t3=t3-ajbk/albd
      if(j.eq.im) t3=t3-aibk/albd
      t4=(w-raa/dam)/am2
      t5=4.d0*w/abm-2.d0*rbb/(abm*dbm)-2.d0*raa/(abm*dam)
      if(i.eq.j.and.ik.eq.im) t5=t5+ctta/(ala*alb)
      if(i.eq.j) t5=t5+rkm/alae-2.d0*bb*ctta/(ala*bm)
      if(ik.eq.im) t5=t5+rij/albe-2.d0*aa*ctta/(alb*am)
      if(i.eq.ik) t5=t5+ajbm/albe
      if(j.eq.ik) t5=t5+aibm/albe
      if(i.eq.im) t5=t5+ajbk/albe
      if(j.eq.im) t5=t5+aibk/albe
      t6=(w-rbb/dbm)/bm2
      t7=(-2.d0*w+raa/dam)/(abm*am)
      if(ik.eq.im) t7=t7+aa*ctta/(alb*am2)
      t8=(-2.d0*w+rbb/dbm)/(abm*bm)
      if(i.eq.j) t8=t8+bb*ctta/(ala*bm2)
      t9=w/(abm*abm)
      go to 140
  130 pl=1.5d0*ctta*ctta-0.5d0
      aibj=a(i)*b(j)
      ajbi=a(j)*b(i)
      akbm=a(ik)*b(im)
      ambk=a(im)*b(ik)
      w=3.d0*ctta
      aibm=a(i)*b(im)*w
      ajbm=a(j)*b(im)*w
      aibk=a(i)*b(ik)*w
      ajbk=a(j)*b(ik)*w
      rkm=akbm+ambk
      rij=aibj+ajbi
      dam=am*gb
      dbm=bm*ga
      raa=w*aa*rkm/dam
      rbb=w*bb*rij/dbm
      aiam=3.d0*a(i)*a(im)/dam
      ajam=3.d0*a(j)*a(im)/dam
      aiak=3.d0*a(i)*a(ik)/dam
      ajak=3.d0*a(j)*a(ik)/dam
      bibm=3.d0*b(i)*b(im)/dbm
      bjbm=3.d0*b(j)*b(im)/dbm
      bibk=3.d0*b(i)*b(ik)/dbm
      bjbk=3.d0*b(j)*b(ik)/dbm
      ama=am*ga
      bmb=bm*gb
      pab=ama*bmb
      e=3.d0*rij*rkm/pab
      akam=3.d0*a(ik)*a(im)/dam
      bibj=3.d0*b(i)*b(j)/dbm
      gab=ga*gb
      t1=aabb*(pl+2.d0/ama+2.d0/bmb-8.d0/pab)+aa*akam/dam+bb*bibj/dbm-rb
     1b-raa+e
      if(i.eq.j) t1=t1+bb*(4.d0/bmb-1.d0)/(ala*ga)
      if(ik.eq.im) t1=t1+aa*(4.d0/ama-1.d0)/(alb*gb)
      if(i.eq.j.and.ik.eq.im) t1=t1-2.d0/(ala*alb*gab)
      if(i.eq.ik) t1=t1+(ajbm-ajam-bjbm)/(r*gab)
      if(j.eq.ik) t1=t1+(aibm-aiam-bibm)/(r*gab)
      if(i.eq.im) t1=t1+(ajbk-ajak-bjbk)/(r*gab)
      if(j.eq.im) t1=t1+(aibk-aiak-bibk)/(r*gab)
      if(i.eq.ik.and.j.eq.im) t1=t1+3.d0/(r*r*gab*gab)
      if(i.eq.im.and.j.eq.ik) t1=t1+3.d0/(r*r*gab*gab)
      t2=2.d0*aabb*(-pl-1.d0/ama-2.d0/bmb+4.d0/pab)/am
      t2=t2-2.d0*aa*akam/(am*dam)+w*bb*rij/(ama*bm)+2.d0*w*aa*rkm/(am*
     1dam)-e/am
      if(i.eq.j) t2=t2+(bb*(pl+2.d0/bmb)+(akam-w*rkm)/dam)/ala
      if(ik.eq.im) t2=t2+2.d0*aa*(1.d0-2.d0/ama)/(alb*dam)
      if(i.eq.j.and.ik.eq.im) t2=t2-1.d0/(ala*alb*gb)
      if(i.eq.ik) t2=t2+(ajam-ajbm)/(r*am*gab)
      if(j.eq.ik) t2=t2+(aiam-aibm)/(r*am*gab)
      if(i.eq.im) t2=t2+(ajak-ajbk)/(r*am*gab)
      if(j.eq.im) t2=t2+(aiak-aibk)/(r*am*gab)
      t3=2.d0*aabb*(-pl-1.d0/bmb-2.d0/ama+4.d0/pab)/bm
      t3=t3-2.d0*bb*bibj/(dbm*bm)+w*aa*rkm/(bmb*am)+2.d0*w*bb*rij/(bm*
     1dbm)-e/bm
      if(i.eq.j) t3=t3+2.d0*bb*(1.d0-2.d0/bmb)/(ala*dbm)
      if(ik.eq.im) t3=t3+(aa*(pl+2.d0/ama)+(bibj-w*rij)/dbm)/alb
      if(i.eq.j.and.ik.eq.im) t3=t3-1.d0/(ala*alb*ga)
      if(i.eq.ik) t3=t3+(bjbm-ajbm)/(r*bm*gab)
      if(j.eq.ik) t3=t3+(bibm-aibm)/(r*bm*gab)
      if(i.eq.im) t3=t3+(bjbk-ajbk)/(r*bm*gab)
      if(j.eq.im) t3=t3+(bibk-aibk)/(r*bm*gab)
      t4=aabb*(pl+2.d0/bmb)/am2+aa*(akam-w*rkm)/(am2*dam)
      if(ik.eq.im) t4=t4-aa/(alb*am*dam)
      t5=4.d0*aabb*(pl+1.d0/ama+1.d0/bmb-2.d0/pab)/abm
      t5=t5-(2.d0*(raa+rbb)-e)/abm
      if(i.eq.j) t5=t5+(w*rkm/dam-2.d0*bb*(pl+1.d0/bmb))/(ala*bm)
      if(ik.eq.im) t5=t5+(w*rij/dbm-2.d0*aa*(pl+1.d0/ama))/(alb*am)
      if(i.eq.j.and.ik.eq.im) t5=t5+pl/(ala*alb)
      if(i.eq.ik) t5=t5+ajbm/(r*abm*gab)
      if(j.eq.ik) t5=t5+aibm/(r*abm*gab)
      if(i.eq.im) t5=t5+ajbk/(r*abm*gab)
      if(j.eq.im) t5=t5+aibk/(r*abm*gab)
      t6=aabb*(pl+2.d0/ama)/bm2+bb*(bibj-w*rij)/(bm2*dbm)
      if(i.eq.j) t6=t6-bb/(ala*dbm*bm)
      t7=(raa-2.d0*aabb*(pl+1.d0/bmb))/(am*abm)
      if(ik.eq.im) t7=t7+aa*pl/(alb*am2)
      t8=(rbb-2.d0*aabb*(pl+1.d0/ama))/(bm*abm)
      if(i.eq.j) t8=t8+bb*pl/(ala*bm2)
      t9=aabb*pl/(abm*abm)
  140 continue
      dpld=t1*f1+t2*f2+t3*f3+t4*f4+t5*f5+t6*f6+t7*f7+t8*f8+t9*f9
      psdint(ij)=dpld
  150 continue
      go to 1000
c......
  200 if(l-1) 210,1000,220
  210 call intga(gb,11)
      bm2=bm*bm
      alb=2.d0*alfb*bm
      ij=0
      do 215 m=1,3
      i=id(m)
      j=jd(m)
      if(iandj) maxk=m
      do 215 mk=1,maxk
      ik=id(mk)
      im=jd(mk)
      ij=ij+1
      bb=b(ik)*b(im)
      t1=bb/r
      t2=-2.d0*bb/(bm*r)
      if(ik.eq.im) t2=t2+1.d0/(alb*r)
      t3=bb/(r*bm2)
      dpld=(t1*f1+t2*f2+t3*f3)*0.33333333333333333d0
      psdint(ij)=dpld
  215 continue
      go to 1000
c....
  220 call intga(gb,12)
      bm2=bm*bm
      alb=2.d0*alfb*bm
      ij=0
      do 225 m=1,6
      i=id(m)
      j=jd(m)
      if(iandj) maxk=m
      do 225 mk=1,maxk
      ik=id(mk)
      im=jd(mk)
      ij=ij+1
      bb=b(ik)*b(im)
      w=3.d0*bb*b(i)*b(j)/bm2
      bmb=bm*gb
      rbb=r*bmb
      t1=w
      if(i.eq.j) t1=t1+bb*(4.d0/bmb-1.d0)/r
      if(i.eq.ik) t1=t1-3.d0*b(j)*b(im)/rbb
      if(i.eq.im) t1=t1-3.d0*b(j)*b(ik)/rbb
      if(j.eq.ik) t1=t1-3.d0*b(i)*b(im)/rbb
      if(j.eq.im) t1=t1-3.d0*b(i)*b(ik)/rbb
      if(i.eq.j.and.ik.eq.im) t1=t1-1.d0/(rbb*alfb)
      if(i.eq.ik.and.j.eq.im) t1=t1+3.d0/(r*r*gb*gb)
      if(j.eq.ik.and.i.eq.im) t1=t1+3.d0/(r*r*gb*gb)
      t2=-2.d0*w/bm
      if(i.eq.j) t2=t2+2.d0*bb*(1.d0-2.d0/bmb)/(r*bm)
      if(ik.eq.im) t2=t2+3.d0*b(i)*b(j)/(bm2*alb)
      if(i.eq.ik) t2=t2+3.d0*b(j)*b(im)/(rbb*bm)
      if(i.eq.im) t2=t2+3.d0*b(j)*b(ik)/(rbb*bm)
      if(j.eq.ik) t2=t2+3.d0*b(i)*b(im)/(rbb*bm)
      if(j.eq.im) t2=t2+3.d0*b(i)*b(ik)/(rbb*bm)
      if(i.eq.j.and.ik.eq.im) t2=t2-1.d0/(r*alb)
      t3=w/bm2
      if(i.eq.j) t3=t3-bb/(r*bm2)
      dpld=(t1*f1+t2*f2+t3*f3)/15.0d0
      psdint(ij)=dpld
  225 continue
      go to 1000
c......
  300 if(l-1) 310,1000,320
  310 call intga(ga,11)
      ala=2.d0*alfa*am
      am2=am*am
      ij=0
      do 315 m=1,6
      i=id(m)
      j=jd(m)
      if(iandj)maxk=m
      do 315 mk=1,maxk
      ik=id(mk)
      im=jd(mk)
      ij=ij+1
      if(ik.ne.im) go to 315
      aa=a(i)*a(j)
      t1=aa/r
      t2=-2.d0*aa/(am*r)
      if(i.eq.j) t2=t2+1.d0/(ala*r)
      t3=aa/(r*am2)
      dpld=(t1*f1+t2*f2+t3*f3)*0.33333333333333333d0
      psdint(ij)=dpld
  315 continue
      go to 1000
  320 call intga(ga,12)
      ala=2.d0*alfa*am
      am2=am*am
      ij=0
      do 325 m=1,6
      i=id(m)
      j=jd(m)
      if(iandj)maxk=m
      do 325 mk=1,maxk
      ik=id(mk)
      im=jd(mk)
      ij=ij+1
      aa=a(i)*a(j)
      w=3.d0*a(ik)*a(im)*aa/am2
      ama=am*ga
      raa=r*ama
      t1=w
      if(ik.eq.im) t1=t1+aa*(4.d0/ama-1.d0)/r
      if(i.eq.ik) t1=t1-3.d0*a(j)*a(im)/raa
      if(i.eq.im) t1=t1-3.d0*a(j)*a(ik)/raa
      if(j.eq.ik) t1=t1-3.d0*a(i)*a(im)/raa
      if(j.eq.im) t1=t1-3.d0*a(i)*a(ik)/raa
      if(i.eq.j.and.ik.eq.im) t1=t1-1.d0/(r*alfa*ama)
      if(i.eq.ik.and.j.eq.im) t1=t1+3.d0/(r*r*ga*ga)
      if(j.eq.ik.and.i.eq.im) t1=t1+3.d0/(r*r*ga*ga)
      t2=-2.d0*w/am
      if(ik.eq.im) t2=t2+2.d0*aa*(1.d0-2.d0/ama)/(r*am)
      if(i.eq.j) t2=t2+3.d0*a(ik)*a(im)/(ala*am2)
      if(i.eq.ik) t2=t2+3.d0*a(j)*a(im)/(raa*am)
      if(i.eq.im) t2=t2+3.d0*a(j)*a(ik)/(raa*am)
      if(j.eq.ik) t2=t2+3.d0*a(i)*a(im)/(raa*am)
      if(j.eq.im) t2=t2+3.d0*a(i)*a(ik)/(raa*am)
      if(i.eq.j.and.ik.eq.im) t2=t2-1.d0/(r*ala)
      t3=w/am2
      if(ik.eq.im) t3=t3-aa/(r*am2)
      dpld=(t1*f1+t2*f2+t3*f3)/15.0d0
      psdint(ij)=dpld
  325 continue
      go to 1000
c......
  400 if(l-1) 410,1000,420
  410 dpld=fgam(7-n)/(18.0d0*r*r)
      ij=0
      do 415 m=1,3
      if(iandj)maxk=m
      do 415 mk=1,maxk
      ij=ij+1
      if(mk.gt.3) go to 415
      psdint(ij)=dpld
  415 continue
      go to 1000
  420 t1=fgam(7-n)/(450.0d0*r*r)
      ij=0
      do 425 m=1,6
      i=id(m)
      j=jd(m)
      if(iandj) maxk=m
      do 425 mk=1,maxk
      ik=id(mk)
      im=jd(mk)
      ij=ij+1
      dpld=0.d0
      if(i.eq.j.and.ik.eq.im) dpld=dpld-2.0d0*t1
      if(i.eq.ik.and.j.eq.im) dpld=dpld+3.d0*t1
      if(j.eq.ik.and.i.eq.im) dpld=dpld+3.d0*t1
      psdint(ij)=dpld
  425 continue
      go to 1000
c......
  500 go to (600,700,800,900),k
  600 if(cm.lt.1.d-8) go to 670
      s=1.d0/cm
      w=0.5d0/rab
      call intga(gc,22)
      ij=0
      do 605 m=1,6
      i=id(m)
      j=jd(m)
      if(iandj) maxk=m
      do 605 mk=1,maxk
      ik=id(mk)
      im=jd(mk)
      ij=ij+1
      aa=a(i)*a(j)
      bb=b(ik)*b(im)
      aabb=aa*bb
      ccij=c(i)*c(j)
      cckm=c(ik)*c(im)
      caac=c(i)*a(j)+a(i)*c(j)
      cbbc=c(ik)*b(im)+b(ik)*c(im)
      t1=aabb
      t2=-aa*cbbc-bb*caac
      t3=aa*cckm+bb*ccij+caac*cbbc
      t4=-ccij*cbbc-cckm*caac
      t5=ccij*cckm
      if(i.ne.j) go to 610
      t2=t2+w*bb
      t3=t3-w*cbbc
      if(ik.eq.im) t3=t3+w*w
      t4=t4+w*cckm
  610 if(ik.ne.im) go to 620
      t2=t2+w*aa
      t3=t3-w*caac
      t4=t4+w*ccij
  620 if(i.ne.ik) go to 630
      t2=t2+w*a(j)*b(im)
      t3=t3-w*(c(j)*b(im)+a(j)*c(im))
      if(j.eq.im) t3=t3+w*w
      t4=t4+w*c(j)*c(im)
  630 if(j.ne.im) go to 640
      t2=t2+w*a(i)*b(ik)
      t3=t3-w*(a(i)*c(ik)+c(i)*b(ik))
      t4=t4+w*c(i)*c(ik)
  640 if(i.ne.im) go to 650
      t2=t2+w*a(j)*b(ik)
      t3=t3-w*(a(j)*c(ik)+c(j)*b(ik))
      if(j.eq.ik) t3=t3+w*w
      t4=t4+w*c(j)*c(ik)
  650 if(j.ne.ik) go to 660
      t2=t2+w*a(i)*b(im)
      t3=t3-w*(c(i)*b(im)+a(i)*c(im))
      t4=t4+w*c(i)*c(im)
  660 continue
      dpld=(((t5*f5*s+t4*f4)*s+t3*f3)*s+t2*f2)*s+t1*f1
      psdint(ij)=dpld
  605 continue
      go to 1000
c....
  670 dxpt=dexp(-pt)
      t1=fgam(5-n)/(6.d0*r)
      t2=fgam(7-n)/(30.d0*r*r)
      ij=0
      do 675 m=1,6
      i=id(m)
      j=jd(m)
      if(iandj) maxk=m
      do 675 mk=1,maxk
      ik=id(mk)
      im=jd(mk)
      ij=ij+1
      s=0.0d0
      if(i.eq.j) s=s+b(ik)*b(im)
      if(ik.eq.im) s=s+a(i)*a(j)
      if(i.eq.ik) s=s+a(j)*b(im)
      if(i.eq.im) s=s+a(j)*b(ik)
      if(j.eq.ik) s=s+a(i)*b(im)
      if(j.eq.im) s=s+a(i)*b(ik)
      dpld=0.5d0*a(i)*a(j)*b(ik)*b(im)*fgam(3-n)
      dpld=dpld+s*t1
      s=t2
      if((i.eq.j).and.(ik.eq.im)) dpld=dpld+s
      if((i.eq.ik).and.(j.eq.im)) dpld=dpld+s
      if((i.eq.im).and.(j.eq.ik)) dpld=dpld+s
      psdint(ij)=dpld*dxpt
  675 continue
      go to 1000
c......
  700 if(cm.lt.1.d-8) go to 770
      s=1.d0/cm
      w=0.5d0/rab
      call intga(gc,23)
      ij=0
      do 705 m=1,6
      i=id(m)
      j=jd(m)
      if(iandj)maxk=m
      do 705 mk=1,maxk
      ik=id(mk)
      im=jd(mk)
      ij=ij+1
      bb=b(ik)*b(im)
      ccij=c(i)*c(j)
      cckm=c(ik)*c(im)
      cbbc=c(ik)*b(im)+b(ik)*c(im)
      t1=0.0d0
      t2=bb*ccij
      t3=-ccij*cbbc
      t4=ccij*cckm
      if(i.ne.j) go to 710
      t1=t1+w*bb
      t2=t2-w*cbbc
      if(ik.eq.im) t2=t2+w*w
      t3=t3+w*cckm
  710 if(ik.ne.im) go to 720
      t3=t3+w*ccij
  720 if(i.ne.ik) go to 730
      t2=t2-w*c(j)*b(im)
      if(j.eq.im) t2=t2+w*w
      t3=t3+w*c(j)*c(im)
  730 if(j.ne.im) go to 740
      t2=t2-w*c(i)*b(ik)
      t3=t3+w*c(i)*c(ik)
  740 if(i.ne.im) go to 750
      t2=t2-w*c(j)*b(ik)
      if(j.eq.ik) t2=t2+w*w
      t3=t3+w*c(j)*c(ik)
  750 if(j.ne.ik) go to 760
      t2=t2-w*c(i)*b(im)
      t3=t3+w*c(i)*c(im)
  760 continue
      dpld=(((t4*f4*s+t3*f3)*s+t2*f2)*s+t1*f1)*s
      psdint(ij)=dpld
  705 continue
      go to 1000
  770 dxpt=dexp(-pt)
      t1=fgam(7-n)/(30.d0*r*r)
      t2=fgam(5-n)/(6.d0*r)
      ij=0
      do 775 m=1,6
      i=id(m)
      j=jd(m)
      if(iandj) maxk=m
      do 775 mk=1,maxk
      ik=id(mk)
      im=jd(mk)
      ij=ij+1
      dpld=0.d0
      if((i.eq.j).and.(ik.eq.im)) dpld=dpld+1.0d0
      if((i.eq.ik).and.(j.eq.im)) dpld=dpld+1.0d0
      if((i.eq.im).and.(j.eq.ik)) dpld=dpld+1.0d0
      dpld=dpld*t1
      if(i.eq.j) dpld=dpld+b(ik)*b(im)*t2
      dpld=dpld*dxpt
      psdint(ij)=dpld
  775 continue
      go to 1000
c......
  800 if(cm.lt.1.d-8) go to 870
      s=1.d0/cm
      w=0.5d0/rab
      call intga(gc,23)
      ij=0
      do 805 m=1,6
      i=id(m)
      j=jd(m)
      if(iandj) maxk=m
      do 805 mk=1,maxk
      ik=id(mk)
      im=jd(mk)
      ij=ij+1
      ccij=c(i)*c(j)
      cckm=c(ik)*c(im)
      caac=c(i)*a(j)+a(i)*c(j)
      t1=0.0d0
      t2=aa*cckm
      t3=-cckm*caac
      t4=ccij*cckm
      if(i.ne.j) go to 810
      t3=t3+w*cckm
      if(ik.eq.im) t2=t2+w*w
  810 if(ik.ne.im) go to 820
      t1=t1+w*aa
      t2=t2-w*caac
      t3=t3+w*ccij
  820 if(i.ne.ik) go to 830
      t2=t2-w*a(j)*c(im)
      if(j.eq.im) t2=t2+w*w
      t3=t3+w*c(j)*c(im)
  830 if(j.ne.im) go to 840
      t2=t2-w*a(i)*c(ik)
      t3=t3+w*c(i)*c(ik)
  840 if(i.ne.im) go to 850
      t2=t2-w*a(j)*c(ik)
      if(j.eq.ik) t2=t2+w*w
      t3=t3+w*c(j)*c(ik)
  850 if(j.ne.ik) go to 860
      t2=t2-w*a(i)*c(im)
      t3=t3+w*c(i)*c(im)
  860 continue
      dpld=(((t4*f4*s+t3*f3)*s+t2*f2)*s+t1*f1)*s
      psdint(ij)=dpld
  805 continue
      go to 1000
  870 dxpt=dexp(-pt)
      t1=fgam(7-n)/(30.d0*r*r)
      t2=fgam(5-n)/(6.d0*r)
      ij=0
      do 875 m=1,6
      i=id(m)
      j=jd(m)
      if(iandj)maxk=m
      do 875 mk=1,maxk
      ik=id(mk)
      im=jd(mk)
      ij=ij+1
      dpld=0.0d0
      if((i.eq.j).and.(ik.eq.im)) dpld=dpld+1.0d0
      if((i.eq.ik).and.(j.eq.im)) dpld=dpld+1.0d0
      if((i.eq.im).and.(j.eq.ik)) dpld=dpld+1.0d0
      dpld=dpld*t1
      if(ik.eq.im) dpld=dpld+a(i)*a(j)*t2
      dpld=dpld*dxpt
      psdint(ij)=dpld
  875 continue
      go to 1000
c......
  900 t1=fgam(7-n)/(30.d0*r*r)
      ij=0
      do 905 m=1,6
      i=id(m)
      j=jd(m)
      if(iandj) maxk=m
      do 905 mk=1,maxk
      ik=id(mk)
      im=jd(mk)
      ij=ij+1
      dpld=0.d0
      if((i.eq.j).and.(ik.eq.im)) dpld=dpld+1.0d0
      if((i.eq.ik).and.(j.eq.im)) dpld=dpld+1.0d0
      if((i.eq.im).and.(j.eq.ik)) dpld=dpld+1.0d0
      dpld=dpld*t1
      psdint(ij)=dpld
  905 continue
c......
 1000 continue
      return
      end
