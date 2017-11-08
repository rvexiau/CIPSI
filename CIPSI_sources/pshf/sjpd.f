      subroutine sjpd(nshela,numa,katoa,ca,kstara,kna,ktypa,kloca,
     *                exa,csa,cpa,cda,cfa,cga,
     *                nshelb,numb,katob,cb,kstarb,knb,ktypb,klocb,
     *                exb,csb,cpb,cdb,cfb,cgb,
     *                itol,icut,r)
      implicit real*8 (a-h,o-z)
      logical iandj,norm,double
      common/stv/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      dimension katoa(*),ca(3,*),kstara(*),kna(*),ktypa(*),
     * kloca(*),exa(*),csa(*),cpa(*),cda(*),cfa(*),
     *cga(*)
      dimension katob(*),cb(3,*),kstarb(*),knb(*),ktypb(*),
     * klocb(*),exb(*),csb(*),cpb(*),cdb(*),cfb(*),
     *cgb(*)
      dimension r(*)
      dimension s(225),dij(225),
     1 xin(125),yin(125),zin(125)
      dimension ix(35),iy(35),iz(35),jx(35),jy(35),jz(35)
      dimension ijx(225),ijy(225),ijz(225)
      dimension kmin(10),kmax(10),kmaz(10)
      data pi212 /1.1283791670955d+00/
      data sqrt3/1.73205080756888d+00/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/
c
      data jx / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
     1          3, 0, 0, 2, 2, 1, 0, 1, 0, 1,
     2          4, 0, 0, 3, 3, 1, 0, 1, 0, 2,
     3          2, 0, 2, 1, 1/
      data ix / 1, 6, 1, 1,11, 1, 1, 6, 6, 1,
     1         16, 1, 1,11,11, 6, 1, 6, 1, 6,
     2         21, 1, 1,16,16, 6, 1, 6, 1,11,
     3         11, 1,11, 6, 6/
      data jy / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
     1          0, 3, 0, 1, 0, 2, 2, 0, 1, 1,
     2          0, 4, 0, 1, 0, 3, 3, 0, 1, 2,
     3          0, 2, 1, 2, 1/
      data iy / 1, 1, 6, 1, 1,11, 1, 6, 1, 6,
     1          1,16, 1, 6, 1,11,11, 1, 6, 6,
     2          1,21, 1, 6, 1,16,16, 1, 6,11,
     3          1,11, 6,11, 6/
      data jz / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
     1          0, 0, 3, 0, 1, 0, 1, 2, 2, 1,
     2          0, 0, 4, 0, 1, 0, 1, 3, 3, 0,
     3          2, 2, 1, 1, 2/
      data iz / 1, 1, 1, 6, 1, 1,11, 1, 6, 6,
     1          1, 1,16, 1, 6, 1, 6,11,11, 6,
     2          1, 1,21, 1, 6, 1, 6,16,16, 1,
     3         11,11, 6, 6,11/
      data kmin/1,2,5,11,21,1,1,5,11,21/
      data kmax/1,4,10,20,36,1,4,10,20,36/
      data kmaz/1,4,9,17,28,1,4,10,20,36/
      tol=2.30258d+00*itol
      cutoff=10.d0**(-icut)
c
      norm=.true.
c     ----- ishell
c
      do 9000 ii=1,nshela
      i=katoa(ii)
      xi=ca(1,i)
      yi=ca(2,i)
      zi=ca(3,i)
      i1=kstara(ii)
      i2=i1+kna(ii)-1
      lit=ktypa(ii)
      lis=lit
      if(lis.gt.5) lis=lis-5
      mini=kmin(lit)
      maxi=kmax(lit)
      mazi=kmaz(lit)
      loci=kloca(ii)-mini
c
c     ----- jshell
c
      do 8000 jj=1,nshelb
      j=katob(jj)
      xj=cb(1,j)
      yj=cb(2,j)
      zj=cb(3,j)
      j1=kstarb(jj)
      j2=j1+knb(jj)-1
      ljt=ktypb(jj)
      ljs=ljt
      if(ljs.gt.5) ljs=ljs-5
      minj=kmin(ljt)
      maxj=kmax(ljt)
      mazj=kmaz(ljt)
      locj=klocb(jj)-minj
      nroots=(lis+ljs)/2
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      iandj=.false.
c
c     ----- prepare indices for pairs of (i,j) functions
c
      ij=0
      max=maxj
      do 50 i=mini,maxi
      nx=ix(i)
      ny=iy(i)
      nz=iz(i)
      if(iandj) max=i
      do 50 j=minj,max
      ij=ij+1
      ijx(ij)=nx+jx(j)
      ijy(ij)=ny+jy(j)
      ijz(ij)=nz+jz(j)
   50 continue
      do 60 i=1,ij
      s(i)=0.0d+00
   60 continue
c
c     ----- i primitive
c
      jgmax=j2
      do 7000 ig=i1,i2
      ai=exa(ig)
      arri=ai*rr
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
      csi=csa(ig)
      cpi=cpa(ig)
      cdi=cda(ig)
      cfi=cfa(ig)
      cgi=cga(ig)
c
c
c     ----- j primtive
c
      if(iandj) jgmax=ig
      do 6000 jg=j1,jgmax
      aj=exb(jg)
      aa=ai+aj
      dum=aj*arri/aa
      if(dum.gt.tol) go to 6000
      fac=dexp(-dum)
      csj=csb(jg)
      cpj=cpb(jg)
      cdj=cdb(jg)
      cfj=cfb(jg)
      cgj=cgb(jg)
      ax=(axi+aj*xj)/aa
      ay=(ayi+aj*yj)/aa
      az=(azi+aj*zj)/aa
c
c     ----- density factor
c
      double=iandj.and.ig.ne.jg
      max=maxj
      nn=0
      do 310 i=mini,maxi
      go to ( 70, 80,180,180, 90,180,180,100,180,180,
     1       110,180,180,120,180,180,180,180,180,130,
     2       140,180,180,150,180,180,180,180,180,160,
     3       180,180,170,180,180),i
   70 dum1=csi*fac
      go to 180
   80 dum1=cpi*fac
      go to 180
   90 dum1=cdi*fac
      go to 180
  100 if(norm) dum1=dum1*sqrt3
      go to 180
  110 dum1=cfi*fac
      go to 180
  120 if(norm) dum1=dum1*sqrt5
      go to 180
  130 if(norm) dum1=dum1*sqrt3
      go to 180
  140 dum1=cgi*fac
      go to 180
  150 if(norm) dum1=dum1*sqrt7
      go to 180
  160 if(norm) dum1=dum1*sqrt5/sqrt3
      go to 180
  170 if(norm) dum1=dum1*sqrt3
  180 if(iandj) max=i
      do 310 j=minj,max
      go to (190,200,300,300,210,300,300,220,300,300,
     1       230,300,300,240,300,300,300,300,300,250,
     2       260,300,300,270,300,300,300,300,300,280,
     3       300,300,290,300,300),j
  190 dum2=dum1*csj
      if(.not.double) go to 300
      if(i.gt.1) go to 195
      dum2=dum2+dum2
      go to 300
  195 dum2=dum2+csi*cpj*fac
      go to 300
  200 dum2=dum1*cpj
      if(double) dum2=dum2+dum2
      go to 300
  210 dum2=dum1*cdj
      if(double) dum2=dum2+dum2
      go to 300
  220 if(norm) dum2=dum2*sqrt3
      go to 300
  230 dum2=dum1*cfj
      if(double) dum2=dum2+dum2
      go to 300
  240 if(norm) dum2=dum2*sqrt5
      go to 300
  250 if(norm) dum2=dum2*sqrt3
      go to 300
  260 dum2=dum1*cgj
      if(double) dum2=dum2+dum2
      go to 300
  270 if(norm) dum2=dum2*sqrt7
      go to 300
  280 if(norm) dum2=dum2*sqrt5/sqrt3
      go to 300
  290 if(norm) dum2=dum2*sqrt3
  300 nn=nn+1
  310 dij(nn)=dum2
c
c     ----- overlap
c
      t=dsqrt(aa)
      x0=ax
      y0=ay
      z0=az
      in=-5
      do 340 i=1,lis
      in=in+5
      ni=i
      do 340 j=1,ljs
      jn=in+j
      nj=j
      call stvint
      xin(jn)=xint/t
      yin(jn)=yint/t
      zin(jn)=zint/t
  340 continue
      do 350 i=1,ij
      nx=ijx(i)
      ny=ijy(i)
      nz=ijz(i)
      yz=yin(ny)*zin(nz)
      dum  = yz*xin(nx)
      s(i)=s(i)+dij(i)*dum
  350 continue
c
 6000 continue
 7000 continue
      if((lit.ge.3.and.lit.le.5).or.(ljt.ge.3.and.ljt.le.5))
     1 call comgrp(lis,ljs,iandj,s)
c
c     ----- set up overlap matrice
c
      max=mazj
      nn=0
      do 7500 i=mini,mazi
      li=loci+i
      in=(li-1)*numb
      if(iandj) max=i
      do 7500 j=minj,max
      lj=locj+j
      jn=lj+in
      nn=nn+1
      if(dabs(s(nn)).lt.cutoff) s(nn)=0.d0
      r(jn)=s(nn)
 7500 continue
 8000 continue
 9000 continue
      return
      end
