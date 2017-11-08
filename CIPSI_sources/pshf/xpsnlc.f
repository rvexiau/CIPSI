      subroutine xpsnlc
      implicit real*8(a-h,o-z)
      include 'pshf.prm'
      logical*1 iandj,old_psnl
      character*80 apseud,anam,amolcas
      character*8 apsd
      common/psnloc/apseud(dc),iatno(dc)
c
c
      common/symtry/ invt(48),iso(ds,12),ptr(3,144),dtr(6,288),
     1 ftr(10,480),nt,ict(dc,48)
      common/output/nprint
      common/eigen2/ f(doas),v(7*nopnl,doa)
      common/hsym/ttt(35,35),mini,maxi,lit,minj,maxj,ljt,ntr
      common/iofile/ ir,iw,ip,is
      common/infoa/nat,ich,mul,numscf,nnp,ne,na,nb,zan(dc),c(3,dc)
      common/nshel/ ex(dgp),cs(dgp),cp(dgp),cd(dgp),cf(dgp),cg(dgp),
     1kstart(ds),katom(ds),ktype(ds),kng(ds),kloc(ds),kmin(ds),kmax(ds),
     2nshell,kmaz(ds)
c
      common/isopac/indin(48),indout(12)
      common/recomb/s(225,nopnl),itpmax
      common/stv/ xint,yint,zint,ta,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      dimension m1(doa),h(doas),dij(225),g(doa)
      dimension hh(doas)
      dimension xin(125),yin(125),zin(125),ijx(225),ijy(225),ijz(225),
     1 ix(35),iy(35),iz(35),jx(35),jy(35),jz(35)
      dimension apsl(nopnl),cpsl(nopnl),cpslo(nopnl*nopnl)
      data jx/0,1,0,0,2,0,0,1,1,0,
     1        3,0,0,2,2,1,0,1,0,1,
     2        4,0,0,3,3,1,0,1,0,2,
     3        2,0,2,1,1/
      data ix/1,6,1,1,11,1,1,6,6,1,
     1       16,1,1,11,11,6,1,6,1,6,
     2       21,1,1,16,16,6,1,6,1,11,
     3       11,1,11,6,6/
      data jy/0,0,1,0,0,2,0,1,0,1,
     1        0,3,0,1,0,2,2,0,1,1,
     2        0,4,0,1,0,3,3,0,1,2,
     3        0,2,1,2,1/
      data iy/1,1,6,1,1,11,1,6,1,6,
     1        1,16,1,6,1,11,11,1,6,6,
     2        1,21,1,6,1,16,16,1,1,11,
     3        1,11,6,11,6/
      data jz/0,0,0,1,0,0,2,0,1,1,
     1        0,0,3,0,1,0,1,2,2,1,
     2        0,0,4,0,1,0,1,3,3,1,
     3        2,2,1,1,2/
      data iz/1,1,1,6,1,1,11,1,6,6,
     1        1,1,16,1,6,1,6,11,11,6,
     2        1,1,21,1,6,1,6,16,16,1,
     3        11,11,6,6,11/
      data sqrt3/1.73205080756888d0/,pi32/5.5683279968317d0/
      data sqrt5/2.23606797749979d0/,sqrt7/2.64575131106459d0/
c
      itol=15
      tol=2.30258d0*itol
      ntwd=(nt+3)/4
c     boucle sur les atomes
c
      nxscf=(numscf*(numscf+1))/2
      do 9000 ic=1,nat
c
      anam=apseud(ic)
      if(anam.eq.'') go to 9000
      if(anam.eq.'pssl'.or.anam.eq.'PSSL') go to 9000
c
c     skip if there is an equivalent er with greater index
c

      do 20 it=1,nt
      icd=ict(ic,it)
      if(icd.gt.ic) go to 9000
   20 continue
c
      do 10 i=1,nxscf
   10 f(i)=0.d0
      old_psnl=.true.
      rewind 20
 2105 read(20,end=2115) apsd,npmax
      if(anam(1:8).eq.apsd(1:8)) go to 2150
      do 2110 j=1,npmax
 2110 read(20)
      go to 2105
 2115 continue
      old_psnl=.false.
      rewind 23
 2125 read(23) amolcas,npmax
      write(6,*) amolcas,npmax
      
      npmax=npmax+1
      if(anam(1:12).eq.amolcas(1:12)) go to 2150
      do 2130 j=1,npmax
 2130 read(23)
      go to 2125
 2150 continue
      write(6,*) 'pseudo trouve'
      do 8900 np=1,npmax
      lkt=np
      if(old_psnl)then
      read(20) itpmax,itport,
     *  (apsl(i),i=1,20),
     * (cpslo(i),i=1,400),
     * ( cpsl(i),i=1,20)
      else
      read(23) lsym,itpmax,itport,
     * (cpsl(i),i=1,itport),
     * (apsl(i),i=1,itpmax),
     * (cpslo(i),i=1,itport*itpmax)
      write(6,*) lsym,itpmax,itport
      end if
      if(itpmax.eq.0) go to 8900
c
c
c     symetrie du pseudo non local
      go to (110,120,130,140,150),np
  110 mink=1
      maxk=1
      mazk=1
      go to 200
  120 mink=2
      maxk=4
      mazk=4
      go to 200
  130 mink=5
      maxk=10
      mazk=9
      go to 200
  140 mink=11
      maxk=20
      mazk=17
      go to 200
  150 mink=21
      maxk=35
      mazk=29
  200 continue
      xi=c(1,ic)
      yi=c(2,ic)
      zi=c(3,ic)
      kdif=mazk-mink+1
      do 40 ig=1,itpmax
      ax=apsl(ig)
      ax=ax+ax
      go to (41,42,43,44,45),np
   41 cnorm=pi32/(sqrt(ax)*ax)
      go to 46
   42 cnorm=pi32*0.5d0/(sqrt(ax)*ax*ax)
      go to 46
   43 cnorm=pi32*0.75d0/(sqrt(ax)*ax**3)
      go to 46
   44 cnorm=pi32*1.875d0/(sqrt(ax)*ax**4)
      go to 46
   45 cnorm=pi32*6.5625d0/(sqrt(ax)*ax**5)
   46 cnorm=1.d0/sqrt(cnorm)
      do 40 iip=1,itport
      ipp=(iip-1)*itpmax+ig
      cpslo(ipp)=cpslo(ipp)*cnorm
   40 continue
c
      do 8000 jj=1,nshell
      j=katom(jj)
      xj=c(1,j)
      yj=c(2,j)
      zj=c(3,j)
      j1=kstart(jj)
      j2=j1+kng(jj)-1
      ljt=ktype(jj)
      minj=kmin(jj)
      maxj=kmax(jj)
      mazj=kmaz(jj)
      locj=kloc(jj)-minj
      jdif=maxj-minj+1
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      ij=0
      do 300 i=mink,maxk
      nx=ix(i)
      ny=iy(i)
      nz=iz(i)
      do 300 j=minj,maxj
      ij=ij+1
      ijx(ij)=nx+jx(j)
      ijy(ij)=ny+jy(j)
      ijz(ij)=nz+jz(j)
      do 300 ig=1,itpmax
      s(ij,ig)=0.d0
  300 continue
c
c
      do 7000 ig=1,itpmax
      ai=apsl(ig)
      arri=ai*rr
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
c
      do 6000 jg=j1,j2
      aj=ex(jg)
      aa=ai+aj
      dum=aj*arri/aa
      if(dum.gt.tol) go to 6000
      fac=exp(-dum)
      ax=(axi+aj*xj)/aa
      ay=(ayi+aj*yj)/aa
      az=(azi+aj*zj)/aa
c
      nn=0
      dum=fac
      do 400 i=mink,maxk
c
      if(i.eq.8) dum=dum*sqrt3
      if(i.eq.14) dum=dum*sqrt5
      if(i.eq.20) dum=dum*sqrt3
      if(i.eq.24) dum=dum*sqrt7
      if(i.eq.30) dum=dum*sqrt5/sqrt3
      if(i.eq.33) dum=dum*sqrt3
      do 390 j=minj,maxj
      go to (340,350,380,380,360,380,380,370,380,380,
     1       371,380,380,372,380,380,380,380,380,373,
     2       374,380,380,375,380,380,380,380,380,376,
     3       380,380,377,380,380),j
  340 dum1=cs(jg)
      dum2=dum*dum1
      go to 380
  350 dum1=cp(jg)
      dum2=dum*dum1
      go to 380
  360 dum1=cd(jg)
      dum2=dum*dum1
      go to 380
  370 dum2=dum2*sqrt3
      go to 380
  371 dum1=cf(jg)
      dum2=dum*dum1
      go to 380
  372 dum2=dum2*sqrt5
      go to 380
  373 dum2=dum2*sqrt3
      go to 380
  374 dum1=cg(ig)
      dum2=dum*dum1
      go to 380
  375 dum2=dum2*sqrt7
      go to 380
  376 dum2=dum2*sqrt5/sqrt3
      go to 380
  377 dum2=dum2*sqrt3
  380 nn=nn+1
  390 dij(nn)=dum2
  400 continue
c
      ta=sqrt(aa)
      x0=ax
      y0=ay
      z0=az
      if(lkt.eq.1) go to 240
      x01=x0-xi
      y01=y0-yi
      z01=z0-zi
  240 x02=x0-xj
      y02=y0-yj
      z02=z0-zj
      in=-5
      do 450 i=1,lkt
      in=in+5
      ni=i
      do 450 j=1,ljt
      jn=in+j
      nj=j
      call stvint
      xin(jn)=xint/ta
      yin(jn)=yint/ta
      zin(jn)=zint/ta
  450 continue
      do 460 i=1,ij
      nx=ijx(i)
      ny=ijy(i)
      nz=ijz(i)
      s(i,ig)=s(i,ig)+dij(i)*xin(nx)*yin(ny)*zin(nz)
  460 continue
 6000 continue
 7000 continue
      iandj=.false.
      if((lkt.ge.3.and.lkt.le.5).or.(ljt.ge.3.and.ljt.le.5))
     1 call comps(lkt,ljt,iandj)
      ij=0
      do 480 k=mink,mazk
      do 480 j=minj,mazj
      ij=ij+1
      lj=locj+j
      do 480 io=1,itport
      lk=(io-1)*kdif+k-mink+1
      v(lk,lj)=0.d0
      ipp=(io-1)*itpmax
      do 480 ig=1,itpmax
      cx=cpslo(ipp+ig)
      v(lk,lj)=v(lk,lj)+s(ij,ig)*cx
  480 continue
c
 8000 continue
c
c
      npsxxx=itport*kdif
      itxxx=itport
c
      iijj=0
      do 900 ii=1,numscf
      do 900 jj=1,ii
      iijj=iijj+1
      igs=-kdif-mink+1
      do 800 ij=1,itport
      igs=igs+kdif
      do 800 i=mink,mazk
      i1=igs+i
  800 f(iijj)=f(iijj)+cpsl(ij)*v(i1,ii)*v(i1,jj)
  900 continue
c
c
c
c
c
c
 8900 continue
c
c     symetrization
c
      do 2400 it=1,nt
      if(it.eq.1) go to 2310
      icd=ict(ic,it)
      if(icd.ge.ic) go to 2400
      itm=it-1
      do 2005 it1=1,itm
      icd1=ict(ic,it1)
      if(icd1.eq.icd) go to 2400
 2005 continue
      ntr=invt(it)
c
      do 2010 ish=1,nshell
      do 2012 itr=1,ntwd
 2012 indout(itr)=iso(ish,itr)
      call isoout(nt)
 2010 m1(ish)=indin(it)
c
      do 2300 ii=1,nshell
      lit=ktype(ii)
      mini=kmin(ii)
      maxi=kmax(ii)
      mazi=kmaz(ii)
      loci=kloc(ii)-mini
      id=m1(ii)
      locid=kloc(id)-mini
      do 2300 jj=1,ii
      iandj=ii.eq.jj
      ljt=ktype(jj)
      minj=kmin(jj)
      maxj=kmax(jj)
      mazj=kmaz(jj)
      locj=kloc(jj)-minj
      jd=m1(jj)
      locjd=kloc(jd)-minj
      max=mazj
      do 2200 i=mini,mazi
      lci=loci+i
      llc=(lci*(lci-1))/2
      if(iandj)max=i
      do 2200 j=minj,max
      lcj=locj+j
      ttt(i,j)=f(llc+lcj)
      if(iandj) ttt(j,i)=ttt(i,j)
 2200 continue
      call rhr
c
      do 2250 i=mini,mazi
      lcid=locid+i
c
      if(iandj) max=i
      do 2250 j=minj,max
      lcjd=locjd+j
      if(lcid.lt.lcjd) go to 2210
      llc=(lcid*(lcid-1))/2+lcjd
      go to 2220
 2210 llc=(lcjd*(lcjd-1))/2+lcid
 2220 h(llc)=ttt(i,j)
 2250 continue
 2300 continue
      go to 2320
 2310 continue
      icd=ic
      do 2312 i=1,nxscf
 2312 h(i)=f(i)
 2320 continue
c     nprint=1
      if(nprint.ne.1) go to 3333
      write(iw,1000)icd
      call fout(h,numscf)
 3333 continue
 1000 format(//,' pseudopotential integrals /atom', i3)
 1010 format(10f12.8)
c
      call reada (hh,nxscf,3)
c
      do 2330 i=1,nxscf
 2330 hh(i)=hh(i)+h(i)
      call wrtda(hh,nxscf,3)
 2400 continue
c
 9000 continue
      if(nprint.eq.1) then
      call reada(hh,nxscf,3)
      write(iw,7603)
      call fout(hh,numscf)
 7603 format(/,' one electron integrals (total)')
      end if
      return
c
      end
