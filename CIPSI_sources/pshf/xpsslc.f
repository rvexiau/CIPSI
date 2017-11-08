      subroutine xpsslc
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      logical spps
      character*80 apseud
      logical*1 iandj,norm,out,double
      logical pseud
      common/psnloc/apseud(dc),iatno(dc)
      common /infpot/apot(400),cpot(400),npot(400),nbtyp(200),
     1ipseud(103),nsom(50),pseud
      common/output/nprint,itol,icut,normf,normp
      common/iofile/ir,iw,ip,is
      common/infoa/nat,ich,mul,num,nnp,ne,na,nb,zan(dc),c(3,dc)
      common/nshel/ex(dgp),cs(dgp),cp(dgp),cd(dgp),cf(dgp),cg(dgp),
     1kstart(ds),katom(ds),ktype(ds),kng(ds),kloc(ds),kmin(ds),kmax(ds),
     2nshell,kmaz(ds)
c     common/eigen2/q(doas),r(doas),s(36),dij(36)
      common/eigen2/q(doas),r(doas),s(225),dij(225)
c     common/pmat/ psdint(36),pspt(36),iandj
      common/pmat/ psdint(225),pspt(225),iandj
      common/pdpot/ xi,yi,zi,icent,ni,ai,
     1              xj,yj,zj,jcent,nj,aj,
     2              cx,cy,cz,ic
      data pi212 /1.1283791670955d+00/
      data sqrt3/1.73205080756888d+00/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/
      ia(i)=i*(i-1)/2
      if(.not.pseud) return
c
      call psepin
      write(iw,9999)
 9999 format(//,10x,20(1h*),//,10x,'effective potential(semi-local) inte
     1grals',
     2 //,10x,20(1h*))
      tol=2.30258d+00*itol
      cutoff=10.d0**(-icut)
      out=nprint.eq.1
      norm=normf.ne.1.or.normp.ne.1
      nxscf=num*(num+1)/2
      call reada(q,nxscf,3)
c
      do 10000 ic=1,nat
      write(6,*) 'ic,iatno(ic),ipseud(iatno(ic)',
     *            ic,iatno(ic),ipseud(iatno(ic))
      if(ipseud(iatno(ic)).eq.0)go to 10000
c     ----- ishell
c
      cx=c(1,ic)
      cy=c(2,ic)
      cz=c(3,ic)
      do 9000 ii=1,nshell
      i=katom(ii)
      icent=i
      xi=c(1,i)
      yi=c(2,i)
      zi=c(3,i)
      i1=kstart(ii)
      i2=i1+kng(ii)-1
      lit=ktype(ii)
      lis=lit
      if(lis.gt.5) lis=lis-5
      ni=lis
      mini=kmin(ii)
      maxi=kmax(ii)
      mazi=kmaz(ii)
      loci=kloc(ii)-mini
c
c     ----- jshell
c
      do 8000 jj=1,ii
      j=katom(jj)
      jcent=j
      xj=c(1,j)
      yj=c(2,j)
      zj=c(3,j)
      j1=kstart(jj)
      j2=j1+kng(jj)-1
      ljt=ktype(jj)
      ljs=ljt
      if(ljs.gt.5) ljs=ljs-5
      nj=ljs
      minj=kmin(jj)
      maxj=kmax(jj)
      mazj=kmaz(jj)
      locj=kloc(jj)-minj
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      iandj=ii.eq.jj
c
c     ----- prepare indices for pairs of (i,j) functions
c
      ij=(maxi-mini+1)*(maxj-minj+1)
      if(iandj)ij=(maxi-mini+1)*(maxj-minj+2)/2
      do  i=1,ij
      s(i)=0.0d+00
      end do
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
      dum=pi212/aa
      facinv=aa/(fac*pi212)
      do 380 i=1,ij
  380 dij(i)=dij(i)*dum
c
      jdil=maxj-minj+1
      inin=mini
      inax=maxi
      jnin=minj
      jnax=maxj
c     particular case of sp shells
c
      if(lit-7) 2400,2410,2420
 2400 if(ljt-7) 2490,2430,2490
 2410 if(ljt-7) 2440,2450,2470
 2420 if(ljt-7) 2490,2480,2490
c     s / sp
c     first calculate s /s
 2430 nj =1
      jnin=1
      jnax=1
      iret=1
      go to 2500
c     second calculate s / p
 2435 nj =2
      jnin=2
      jnax=4
      iret=8
      go to 2500
c     sp / s
c     first calculate s / s
 2440 ni =1
      inin=1
      inax=1
      iret=2
      go to 2500
c     second calculate p /s
 2445 ni =2
      inin=2
      inax=4
      iret=8
      go to 2500
c     sp / sp
c first calculate s / s
 2450 ni =1
      inin=1
      inax=1
      nj =1
      jnin=1
      jnax=1
      iret=3
      go to 2500
c     second calculate p / s
 2455 ni =2
      inin=2
      inax=4
      iret=4
      go to 2500
c     third calculate s / p
 2460 nj=2
      jnin=2
      jnax=4
      if(iandj) go to 2465
      ni=1
      inin=1
      inax=1
      iret=5
      go to 2500
c     fourth calculate p / p
 2465 ni =2
      inin=2
      inax=4
      iret=8
      go to 2500
c     sp / d
c     first calculate s / d
 2470 ni =1
      inin=1
      inax=1
      iret=6
      go to 2500
c     second calculate p / d
 2475 ni =2
      inin=2
      inax=4
      iret=8
      go to 2500
c     d / sp
c     first calculate d / s
 2480 nj =1
      jnin=1
      jnax=1
       iret=7
      go to 2500
c     second calculate d/p
 2485 nj =2
      jnin=2
      jnax=4
      iret=8
      go to 2500
c     general case without sp shells
 2490 iret=8
 2500 spps=iandj.and.iret.eq.8
      nprim=(inax-inin+1)*(jnax-jnin+1)
      if(spps) nprim=ia(inax-inin+2)
      call psepot(nprim)
      idif=inax-inin+1
      jdif=jnax-jnin+1
      nax=jnax
      do 450 i=inin,inax
      if(spps) nax=i
      do 450 j=jnin,nax
      nn=(i-mini)  *jdil+j-minj+1
      if(iandj) nn=ia(i-mini+1)+j-minj+1
      if(spps) go to 446
      if(ni.lt.nj) go to 444
      ips=(i-inin)*jdif+(j-jnin+1)
      go to 448
  444 ips=(j-jnin)*idif+(i-inin+1)
      go to 448
  446 ips=ia(i-inin+1)+(j-jnin+1)
c
  448 s(nn)=s(nn)+pspt(ips)*dij(nn)*facinv
  450 continue
      go to (2435,2445,2455,2460,2465,2475,2485,500),iret
  500 continue
 6000 continue
 7000 continue
      if((lit.ge.3.and.lit.le.5).or.(ljt.ge.3.and.ljt.le.5))
     1 call comone(lis,ljs,iandj)
c
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
      r(jn)=s(nn)
 7500 continue
 7600 format(2i5,2e20.12)
 8000 continue
 9000 continue
      if(out) then
      write(iw,7601)ic
      call fout (r,num)
 7601 format(/,' effective potential(semi-local)integrals/atom ',i4)
      end if
c
      do i=1,nxscf
      q(i)=q(i)+r(i)
      end do
10000 continue
      call wrtda(q,nxscf,3)
      write(iw,9400)
 9400 format(' ...... end of effective potential(semi-local)integrals ..
     1....')
      return
      end
