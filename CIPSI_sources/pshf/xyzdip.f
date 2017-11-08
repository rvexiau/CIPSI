      subroutine xyzdip
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      logical*1 iandj,out,norm,double
      common/output/nprint,itol,icut,normf,normp,nopk
      common/iofile/ir,iw,ip,is,iq,ih,iv,if2,if3,if4,ioda(19)
      common/infoa/nat,ich,mul,num,nnx,ne,na,nb,zan(dc),c(3,dc)
      common/nshel/ex(dgp),cs(dgp),cp(dgp),cd(dgp),cf(dgp),cg(dgp),
     1kstart(ds),katom(ds),ktype(ds),kng(ds),kloc(ds),kmin(ds),kmax(ds),
     2nshell,kmaz(ds)
      common/intdip/xint0,yint0,zint0,xintx,yinty,zintz,t,x0,y0,z0,
     1 xi,yi,zi,xj,yj,zj,ni,nj
c     common/dipxyz/xs(5050),ys(5050),zs(5050)
      common/dipxyz/xs(doas),ys(doas),zs(doas)
      common/recomb/v(225,nopnl),itpmax
      dimension sx(225),sy(225),sz(225)
      equivalence (v(1,1),sx(1)),(v(1,2),sy(1)),(v(1,3),sz(1))
      dimension dij(225),xin(50),yin(50),zin(50)
      dimension ix(35),iy(35),iz(35),jx(35),jy(35),jz(35)
      dimension ijx(225),ijy(225),ijz(225)
      data zero /0.0d+00/
      data sqrt3 /1.73205080756888d+00/
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
      tol=2.30258d+00*itol
      out=nprint.eq.1
      itpmax=3
      norm=normf.ne.1.or.normp.ne.1
      
c
c     ----- ishell
c
      do 9000 ii=1,nshell
      i=katom(ii)
      xi=c(1,i)
      yi=c(2,i)
      zi=c(3,i)
      i1=kstart(ii)
      i2=i1+kng(ii)-1
      lit=ktype(ii)
      lis=lit
      if(lis.gt.5) lis=lis-5
      mini=kmin(ii)
      maxi=kmax(ii)
      mazi=kmaz(ii)
      loci=kloc(ii)-mini
c
c     ----- jshell
c
      do 8000 jj=1,ii
      j=katom(jj)
      xj=c(1,j)
      yj=c(2,j)
      zj=c(3,j)
      j1=kstart(jj)
      j2=j1+kng(jj)-1
      ljt=ktype(jj)
      ljs=ljt
      if(ljs.gt.5) ljs=ljs-5
      minj=kmin(jj)
      maxj=kmax(jj)
      mazj=kmaz(jj)
      locj=kloc(jj)-minj
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      iandj=ii.eq.jj
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
      sx(i)=zero
      sy(i)=zero
   60 sz(i)=zero
c
c     ----- i primitive
c
      jgmax=j2
      do 7000 ig=i1,i2
      ai=ex(ig)
      arri=ai*rr
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
      csi=cs(ig)
      cpi=cp(ig)
      cdi=cd(ig)
      cfi=cf(ig)
      cgi=cg(ig)
c
c
c     ----- j primtive
c
      if(iandj) jgmax=ig
      do 6000 jg=j1,jgmax
      aj=ex(jg)
      aa=ai+aj
      dum=aj*arri/aa
      if(dum.gt.tol) go to 6000
      fac=exp(-dum)
      csj=cs(jg)
      cpj=cp(jg)
      cdj=cd(jg)
      cfj=cf(jg)
      cgj=cg(jg)
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
c
c     ----- dipole moment integrals -----
c
      t=sqrt(aa)
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
      call dipint
      xin(jn)=xint0/t
      yin(jn)=yint0/t
      zin(jn)=zint0/t
      xin(jn+25)=xintx/t
      yin(jn+25)=yinty/t
      zin(jn+25)=zintz/t
  340 continue
      do 350 i=1,ij
      nx=ijx(i)
      ny=ijy(i)
      nz=ijz(i)
      sx(i)=sx(i)+dij(i)*(xin(nx+25)*yin(ny  )*zin(nz  ))
      sy(i)=sy(i)+dij(i)*(xin(nx  )*yin(ny+25)*zin(nz  ))
      sz(i)=sz(i)+dij(i)*(xin(nx  )*yin(ny  )*zin(nz+25))
  350 continue
 6000 continue
 7000 continue
c
      if((lit.ge.3.and.lit.le.5).or.(ljt.ge.3.and.ljt.le.5))
     $call comps(lit,ljt,iandj)
c     ----- set up dipole moment matrices -----
c
      max=mazj
      nn=0
      do 7500 i=mini,mazi
      li=loci+i
      in=(li*(li-1))/2
      if(iandj) max=i
      do 7500 j=minj,max
      lj=locj+j
      jn=lj+in
      nn=nn+1
      xs(jn)=sx(nn)
      ys(jn)=sy(nn)
      zs(jn)=sz(nn)
 7500 continue
 8000 continue
 9000 continue
      write(13) (xs(i),i=1,nnx)
      write(13) (ys(i),i=1,nnx)
      write(13) (zs(i),i=1,nnx)
      if(.not.out) go to 9100
      write(iw,9999)
      call matout(xs,num)
      write(iw,9998)
      call matout(ys,num)
      write(iw,9997)
      call matout(zs,num)
 9999 format(/,10x,25(1h-),/,10x,'x-dipole moment integrals',/,
     1 10x,25(1h-))
 9998 format(/,10x,25(1h-),/,10x,'y-dipole moment integrals',/,
     1 10x,25(1h-))
 9997 format(/,10x,25(1h-),/,10x,'z-dipole moment integrals',/,
     1 10x,25(1h-))
 9100 continue
      return
      end
