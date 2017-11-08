      subroutine standv
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      logical*1 iandj,out,norm,double
      common/scfit/aname,maxit,nconv,npunch
      common/times/ti,tx,tim
      common/restar/timlim,irest,nrec,intloc,ist,jst,kst,lst
      common/output/nprint,itol,icut,normf,normp
      common/stvrt/xx,u(9),w(9),nroots
      common/iofile/ir,iw,ip,is
      common/infoa/nat,ich,mul,num,nnp,ne,na,nb,zan(dc),c(3,dc)
      common/nshel/ex(dgp),cs(dgp),cp(dgp),cd(dgp),cf(dgp),cg(dgp),
     1kstart(ds),katom(ds),ktype(ds),kng(ds),kloc(ds),kmin(ds),kmax(ds),
     2nshell,kmaz(ds)
      common/stv/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      common/eigen2/r(doas),q(doas),s(225),g(225),ft(225),dij(225),
     1 xin(125),yin(125),zin(125)
      dimension ix(35),iy(35),iz(35),jx(35),jy(35),jz(35)
      dimension ijx(225),ijy(225),ijz(225)
      dimension  gg(225),ggw(doas)
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
      write(iw,9999)
 9999 format(//,10x,20(1h*),//,10x,'1 electron integrals',
     2 //,10x,20(1h*))
      tol=2.30258d+00*itol
      cutoff=10.d0**(-icut)
      out=nprint.eq.1
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
      nroots=(lis+ljs)/2
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
      if(j.gt.1) go to 10
      ft(ij)=3.0d+00
      go to 50
   10 if(j.gt.4) go to 20
      ft(ij)=5.0d+00
      go to 50
   20 if(j.gt.10) go to 30
      ft(ij)=7.0d+00
      go to 50
   30 if(j.gt.20) go to 40
      ft(ij)=9.0d+00
      go to 50
   40 ft(ij)=11.0d+00
   50 continue
      do 60 i=1,ij
      s(i)=0.0d+00
      gg(i)=0.d+00
   60 g(i)=0.0d+00
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
c     ----- overlap and kinetic energy
c
      t=sqrt(aa)
      t1=-2.0d+00*aj*aj/t
      t2=-0.5d+00/t
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
      nj=j+2
      call stvint
      xin(jn+25)=xint*t1
      yin(jn+25)=yint*t1
      zin(jn+25)=zint*t1
      nj=j-2
      if(nj.gt.0) go to 320
      xint=0.0d+00
      yint=0.0d+00
      zint=0.0d+00
      go to 330
  320 call stvint
  330 n=(j-1)*(j-2)
      dum=dble(n)*t2
      xin(jn+50)=xint*dum
      yin(jn+50)=yint*dum
      zin(jn+50)=zint*dum
  340 continue
      do 350 i=1,ij
      nx=ijx(i)
      ny=ijy(i)
      nz=ijz(i)
      yz=yin(ny)*zin(nz)
      dum  = yz*xin(nx)
      dum1 = (xin(nx+25)+xin(nx+50))*yz
     1      +(yin(ny+25)+yin(ny+50))*xin(nx)*zin(nz)
     2      +(zin(nz+25)+zin(nz+50))*xin(nx)*yin(ny)
      s(i)=s(i)+dij(i)*dum
      g(i)=g(i)+dij(i)*(dum*aj*ft(i)+dum1)
  350 continue
c
c     ..... nuclear attraction
c
      dum=pi212/aa
      do 400 i=1,ij
  400 dij(i)=dij(i)*dum
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
      do 450 ic=1,nat
      znuc=-zan(ic)
      cx=c(1,ic)
      cy=c(2,ic)
      cz=c(3,ic)
      xx=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
      if(nroots.le.3) call wrt123
      if(nroots.eq.4) call wroot4
      if(nroots.eq.5) call wroot5
      mm=0
      do 420 k=1,nroots
      uu=aa*u(k)
      ww=w(k)*znuc
      tt=aa+uu
      t=sqrt(tt)
      x0=(aax+uu*cx)/tt
      y0=(aay+uu*cy)/tt
      z0=(aaz+uu*cz)/tt
      in=-5+mm
      do 410 i=1,lit
      in=in+5
      ni=i
      do 410 j=1,ljt
      jn=in+j
      nj=j
      call stvint
      xin(jn)=xint
      yin(jn)=yint
      zin(jn)=zint*ww
  410 continue
  420 mm=mm+25
      do 440 i=1,ij
      nx=ijx(i)
      ny=ijy(i)
      nz=ijz(i)
      dum=0.0d+00
      mm=0
      do 430 k=1,nroots
      dum=dum+xin(nx+mm)*yin(ny+mm)*zin(nz+mm)
  430 mm=mm+25
      dumdij=dum*dij(i)
      gg(i)=gg(i)+dumdij
  440 g(i)=g(i)+dumdij
  450 continue
 6000 continue
 7000 continue
      if((lit.ge.3.and.lit.le.5).or.(ljt.ge.3.and.ljt.le.5))
     1 call comone(lis,ljs,iandj)
c
c     ----- set up overlap and h-core matrices
c     ----- write h-core matrix out on tape (is)
c
      max=mazj
      nn=0
      do 7500 i=mini,mazi
      li=loci+i
      in=li*(li-1)/2
      if(iandj) max=i
      do 7500 j=minj,max
      lj=locj+j
      jn=lj+in
      nn=nn+1
      if(abs(s(nn)).lt.cutoff) s(nn)=0.d0
      if(abs(g(nn)).lt.cutoff) g(nn)=0.d0
      r(jn)=s(nn)
      q(jn)=g(nn)
      ggw(jn)=gg(nn)
      if(out.and.ii.eq.jj) write(iw,7600) li,lj,s(nn),g(nn)
 7500 continue
 7600 format(2i5,2e20.12)
 8000 continue
 9000 continue
      if(out) then
      write(iw,7601)
      call fout (r,num)
      write(iw,7602)
      call fout(q,num)
      write(iw,7603)
      call fout(ggw,num)
 7601 format(/,' overlap integrals')
 7602 format(/,' kinetic + nuclear integrals')
 7603 format(/,' nuclear integrals')
      end if
      call wrtda(r,nnp,2)
      call wrtda(q,nnp,3)
      write(iw,9400)
 9400 format(' ...... end of one-electron integrals ......')
      call qmat
      irest=5
c
 9500 continue
      call timit(1)
      if(tim.lt.timlim) go to 9700
      write(iw,9600) timlim,nprint,itol,icut,normf,normp,irest,
     1 ist,jst,kst,lst,nrec,intloc
      write(ip,9600) timlim,nprint,itol,icut,normf,normp,irest,
     1 ist,jst,kst,lst,nrec,intloc
 9600 format(f10.0,10i3,i10,i5)
 9700 continue
      return
      end
