      subroutine genral
      implicit real*8 (a-h,o-z)
      logical*1 iandj,kandl,same,out,norm,double
      common/shlinf/ag(20),csa(20),cpa(20),cda(20),cfa(20),cga(20),
     1              bg(20),csb(20),cpb(20),cdb(20),cfb(20),cgb(20),
     1              cg(20),csc(20),cpc(20),cdc(20),cfc(20),cgc(20),
     1              dg(20),csd(20),cpd(20),cdd(20),cfd(20),cgd(20),
     1              xi,yi,zi,xj,yj,zj,rri,xk,yk,zk,xl,yl,zl,rrk,
     1              nga,ngb,ngc,ngd
      common/gout/gout(10000)
      common/ijpair/aa(400),r(400),x1(400),y1(400),z1(400),dd(6400),
     1 ijd(225)
      common/shlnos/qq4,ii,jj,kk,ll,lit,ljt,lkt,llt,
     1 mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     1 loci,locj,lock,locl,ij,kl,ijkl,nij
      common/root/xx,u(9),w(9),nroots
      common/shlt/tol,cutoff,icount,out
      common/misc/maxll,iandj,kandl,same
      common/output/nprint,itol,icut,normf,normp
      common/xyz/xin(5625),yin(5625),zin(5625)
      common/setint/in(9),kn(9),ni,nj,nk,nl,nmax,mmax
     1,bp01,b00,b10,xcp00,xc00,ycp00,yc00,zcp00,zc00,f00
     2,dxij,dyij,dzij,dxkl,dykl,dzkl
      common/dens/dkl(225),dij(225)
      dimension in1(9)
      data sqrt3 /1.73205080756888d+00/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/
      data pi252/34.986836655250d+00/
      if(ijkl.eq.1) go to 6000
      norm=normf.ne.1.or.normp.ne.1
      ni=lit-1
      nj=ljt-1
      nk=lkt-1
      nl=llt-1
      if(ni.ge.5) ni=ni-5
      if(nj.ge.5) nj=nj-5
      if(nk.ge.5) nk=nk-5
      if(nl.ge.5) nl=nl-5
      dxij=xi-xj
      dyij=yi-yj
      dzij=zi-zj
      dxkl=xk-xl
      dykl=yk-yl
      dzkl=zk-zl
      nmax=ni+nj
      mmax=nk+nl
      max=nmax+1
      do 10 i=1,max
      n=i-1
      if(n.le.ni) in1(i)=125*n+1
      if(n.gt.ni) in1(i)=125*ni+25*(n-ni)+1
   10 continue
      max=mmax+1
      do 20 k=1,max
      n=k-1
      if(n.le.nk) kn(k)=5*n
      if(n.gt.nk) kn(k)=5*nk+n-nk
   20 continue
      ijt=5*ni+nj+1
      klt=5*nk+nl+1
c
c     ----- k primitive
c
      lgmax=ngd
      do 5000 kg=1,ngc
      ak=cg(kg)
      brrk=ak*rrk
      akxk=ak*xk
      akyk=ak*yk
      akzk=ak*zk
      csk=csc(kg)*pi252
      cpk=cpc(kg)*pi252
      cdk=cdc(kg)*pi252
      cfk=cfc(kg)*pi252
      cgk=cgc(kg)*pi252
c
c     ----- l primitive
c
      if(kandl) lgmax=kg
      do 4000 lg=1,lgmax
      al=dg(lg)
      b=ak+al
      bbrrk=al*brrk/b
      if(bbrrk.gt.tol) go to 4000
      csl=csd(lg)
      cpl=cpd(lg)
      cdl=cdd(lg)
      cfl=cfd(lg)
      cgl=cgd(lg)
      xb=(akxk+al*xl)/b
      yb=(akyk+al*yl)/b
      zb=(akzk+al*zl)/b
      bxbk=b*(xb-xk)
      bybk=b*(yb-yk)
      bzbk=b*(zb-zk)
      bxbi=b*(xb-xi)
      bybi=b*(yb-yi)
      bzbi=b*(zb-zi)
c
c     ----- density factor
c
      double=kandl.and.kg.ne.lg
      n=0
      max=maxl
      do 270 k=mink,maxk
      go to( 30, 40,140,140, 50,140,140, 60,140,140,
     1       70,140,140, 80,140,140,140,140,140, 90,
     2      100,140,140,110,140,140,140,140,140,120,
     3      140,140,130,140,140),k
   30 dum1=csk/b
      go to 140
   40 dum1=cpk/b
      go to 140
   50 dum1=cdk/b
      go to 140
   60 if(norm) dum1=dum1*sqrt3
      go to 140
   70 dum1=cfk/b
      go to 140
   80 if(norm) dum1=dum1*sqrt5
      go to 140
   90 if(norm) dum1=dum1*sqrt3
      go to 140
  100 dum1=cgk/b
      go to 140
  110 if(norm) dum1=dum1*sqrt7
      go to 140
  120 if(norm) dum1=dum1*sqrt5/sqrt3
      go to 140
  130 if(norm) dum1=dum1*sqrt3
  140 if(kandl) max=k
      do 270 l=minl,max
      go to(150,160,260,260,170,260,260,180,260,260,
     1      190,260,260,200,260,260,260,260,260,210,
     2      220,260,260,230,260,260,260,260,260,240,
     3      260,260,250,260,260),l
  150 dum2=dum1*csl
      if(.not.double) go to 260
      if(k.gt.1) go to 155
      dum2=dum2+dum2
      go to 260
  155 dum2=dum2+csk*cpl/b
      go to 260
  160 dum2=dum1*cpl
      if(double) dum2=dum2+dum2
      go to 260
  170 dum2=dum1*cdl
      if(double) dum2=dum2+dum2
      go to 260
  180 if(norm) dum2=dum2*sqrt3
      go to 260
  190 dum2=dum1*cfl
      if(double) dum2=dum2+dum2
      go to 260
  200 if(norm) dum2=dum2*sqrt5
      go to 260
  210 if(norm) dum2=dum2*sqrt3
      go to 260
  220 dum2=dum1*cgl
      if(double) dum2=dum2+dum2
      go to 260
  230 if(norm) dum2=dum2*sqrt7
      go to 260
  240 if(norm) dum2=dum2*sqrt5/sqrt3
      go to 260
  250 if(norm) dum2=dum2*sqrt3
  260 n=n+1
  270 dkl(n)=dum2
c
c     ----- pair of i,j primitives
c
      nn=0
      do 3000 n=1,nij
      dum=bbrrk+r(n)
      if(dum.gt.tol) go to 3000
      do 280 i=1,ij
  280 dij(i)=dd(ijd(i)+nn)
      a=aa(n)
      ab=a*b
      aandb=a+b
      expe=exp(-dum)/sqrt(aandb)
      rho=ab/aandb
      xa=x1(n)
      ya=y1(n)
      za=z1(n)
      xx=rho*((xa-xb)**2+(ya-yb)**2+(za-zb)**2)
      axak=a*(xa-xk)
      ayak=a*(ya-yk)
      azak=a*(za-zk)
      axai=a*(xa-xi)
      ayai=a*(ya-yi)
      azai=a*(za-zi)
      c1x=bxbk+axak
      c2x=a*bxbk
      c3x=bxbi+axai
      c4x=b*axai
      c1y=bybk+ayak
      c2y=a*bybk
      c3y=bybi+ayai
      c4y=b*ayai
      c1z=bzbk+azak
      c2z=a*bzbk
      c3z=bzbi+azai
      c4z=b*azai
c
c     ----- roots and weights for quadrature
c
      if(nroots.le.3) then
        call rt123
      endif
      if(nroots.eq.4) call root4
      if(nroots.eq.5) call root5
      if(nroots.ge.6) call droot
      mm=0
      max=nmax+1
c
c     compute two-electron  integrals for each root
c
      do 2000 m=1,nroots
      u2=u(m)*rho
      f00=expe*w(m)
      do 1900 i=1,max
 1900 in(i)=in1(i)+mm
      dum=ab+u2*aandb
      dum2=dum+dum
      bp01=(a+u2)/dum2
      b00=u2/dum2
      b10=(b+u2)/dum2
      xcp00=(u2*c1x+c2x)/dum
      xc00 =(u2*c3x+c4x)/dum
      ycp00=(u2*c1y+c2y)/dum
      yc00 =(u2*c3y+c4y)/dum
      zcp00=(u2*c1z+c2z)/dum
      zc00 =(u2*c3z+c4z)/dum
      call xyzint
 2000 mm=mm+625
c
c     ----- form (i,j//k,l) integrals over functions
c
      call forms
 3000 nn=nn+16
 4000 continue
 5000 continue
c
c     ----- study symmetry and write integrals out on tape (is)
c
      call qout
      return
 6000 call s0000
      return
      end
