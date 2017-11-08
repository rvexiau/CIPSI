      subroutine ijprim
      implicit real*8 (a-h,o-z)
      logical*1 iandj,kandl,same,out,norm
      common/output/nprint,itol,icut,normf,normp
      common/shlt/tol,cutoff,icount,out
      common/shlinf/ag(20),csa(20),cpa(20),cda(20),cfa(20),cga(20),
     1              bg(20),csb(20),cpb(20),cdb(20),cfb(20),cgb(20),
     1              cg(20),csc(20),cpc(20),cdc(20),cfc(20),cgc(20),
     1              dg(20),csd(20),cpd(20),cdd(20),cfd(20),cgd(20),
     1              xi,yi,zi,xj,yj,zj,rri,xk,yk,zk,xl,yl,zl,rrk,
     1              nga,ngb,ngc,ngd
      common/misc/maxll,iandj,kandl,same
      common/shlnos/qq4,ii,jj,kk,ll,lit,ljt,lkt,llt,
     1 mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     1 loci,locj,lock,locl,ij,kl,ijkl,nij
      common/ijpair/a(400),r(400),x1(400),y1(400),z1(400),dij(6400),
     1 ijd(225)
      data sqrt3 /1.73205080756888d+00/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/
      norm=normf.ne.1.or.normp.ne.1
      max=maxj
      n=0
      nn=0
      ni=lit
      if(ni.ge.6) ni=ni-5
      do 50 i=mini,maxi
      go to (10,10,20,20,10,20,20,10,20,20,10,20,20,10,20,20,20,20,20,10
     1      ,10,20,20,10,20,20,20,20,20,10,20,20,10,20,20),i
   10 nm=nn
   20 nn=nm
      if(iandj) max=i
      do 50 j=minj,max
      go to (30,30,40,40,30,40,40,30,40,40,30,40,40,30,40,40,40,40,40,30
     1      ,30,40,40,30,40,40,40,40,40,30,40,40,30,40,40),j
   30 nn=nn+1
   40 n=n+1
   50 ijd(n)=nn
c
c     ----- i primitive
c
      nij=0
      jbmax=ngb
      do 390 ia=1,nga
      ai=ag(ia)
      arri=ai*rri
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
      csi=csa(ia)
      cpi=cpa(ia)
      cdi=cda(ia)
      cfi=cfa(ia)
      cgi=cga(ia)
c
c     ----- j primitive
c
      if(iandj) jbmax=ia
      do 380 jb=1,jbmax
      aj=bg(jb)
      aa=ai+aj
      dum=aj*arri/aa
      if(dum.gt.tol) go to 380
      csj=csb(jb)
      cpj=cpb(jb)
      cdj=cdb(jb)
      cfj=cfb(jb)
      cgj=cgb(jb)
      nm=16*nij
      nn=nm
      nij=nij+1
      r(nij)=dum
      a(nij)=aa
      x1(nij)=(axi+aj*xj)/aa
      y1(nij)=(ayi+aj*yj)/aa
      z1(nij)=(azi+aj*zj)/aa
c
c     ----- density factor
c
      do 310 i=mini,maxi
      go to( 60, 70,310,310, 80,310,310, 90,310,310,
     1      100,310,310,110,310,310,310,310,310,120,
     2      130,310,310,140,310,310,310,310,310,150,
     3      310,310,160,310,310),i
   60 dum1=csi/aa
      go to 170
   70 dum1=cpi/aa
      go to 170
   80 dum1=cdi/aa
      go to 170
   90 if(norm) dum1=dum1*sqrt3
      go to 170
  100 dum1=cfi/aa
      go to 170
  110 if(norm) dum1=dum1*sqrt5
      go to 170
  120 if(norm) dum1=dum1*sqrt3
      go to 170
  130 dum1=cgi/aa
      go to 170
  140 if(norm) dum1=dum1*sqrt7
      go to 170
  150 if(norm) dum1=dum1*sqrt5/sqrt3
      go to 170
  160 if(norm) dum1=dum1*sqrt3
  170 if(iandj) max=i
      do 300 j=minj,max
      go to(180,190,300,300,200,300,300,210,300,300,
     1      220,300,300,230,300,300,300,300,300,240,
     2      250,300,300,260,300,300,300,300,300,270,
     3      300,300,280,300,300),j
  180 dum2=dum1*csj
      go to 290
  190 dum2=dum1*cpj
      go to 290
  200 dum2=dum1*cdj
      go to 290
  210 if(norm) dum2=dum2*sqrt3
      go to 290
  220 dum2=dum1*cfj
      go to 290
  230 if(norm) dum2=dum2*sqrt5
      go to 290
  240 if(norm) dum2=dum2*sqrt3
      go to 290
  250 dum2=dum1*cgj
      go to 290
  260 if(norm) dum2=dum2*sqrt7
      go to 290
  270 if(norm) dum2=dum2*sqrt5/sqrt3
      go to 290
  280 if(norm) dum2=dum2*sqrt3
  290 nn=nn+1
      dij(nn)=dum2
  300 continue
  310 continue
      if(.not.iandj) go to 380
      if(ia.eq.jb) go to 380
      go to (370,320,350,340,330),ni
  320 if(mini.eq.2) go to 370
      dij(nm+2)=dij(nm+2)+csi*cpj/aa
      go to 360
  330 dij(nm+10)=dij(nm+10)+dij(nm+10)
      dij(nm+9)=dij(nm+9)+dij(nm+9)
      dij(nm+8)=dij(nm+8)+dij(nm+8)
      dij(nm+7)=dij(nm+7)+dij(nm+7)
  340 dij(nm+6)=dij(nm+6)+dij(nm+6)
      dij(nm+5)=dij(nm+5)+dij(nm+5)
      dij(nm+4)=dij(nm+4)+dij(nm+4)
  350 dij(nm+2)=dij(nm+2)+dij(nm+2)
  360 dij(nm+3)=dij(nm+3)+dij(nm+3)
  370 dij(nm+1)=dij(nm+1)+dij(nm+1)
  380 continue
  390 continue
      return
      end
