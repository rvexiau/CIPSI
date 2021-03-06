      subroutine deter(yspipro)
      implicit real*8(a-h,o-x,z),logical*1(y)
      include 'pshf.prm'
      integer*4 trou,part,nexst,nd
      integer*4 iorb,iwspin,itsym,its,isytr
      common/trsf/nt(8),np(8),mt(8),mp(8)
      common /pr/ ncfpr
      common/spo/ feff(2*doa),fdiag(2*doa),
     * iorb(2*doa),iwspin(2*doa),ispin(2*doa),
     *itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),isytr(2*doa,nsymz)
      common/dinde/trou(8*ndimh),part(8*ndimh),nexst(ndimh),
     * nd(ndimh)
      common/e0cip/ebmp(metz),eben(metz),tau
      common  /cipsi/ norb,mnorb,noca,nocb,
     *igela,igelb,ncf,ncper,nsym,isym,isz,ntrsy,metat,
     *isymat(smax),nelac,nprint,yprt,yion
      common/mod90/noa,ydp
      dimension nts(8*lungo),nps(8*lungo),nexs(lungo),nds(lungo)      
      dimension isydeg(nsymz,nsymz),ind(2*doa),
     * itsyp(2*doa)
      dimension indis(doa,nsymz)
      integer*4 iplus,imoins
      data iplus,imoins/1,-1/
cCc       data iplus,imoins/'+','-'/
      igela=0
      maxco=8
      igelb=0
      nsym=1
      ntyp=0
      ntrsy=0
      read(40) nsym,norb,noc,ntrans,(itsyp(i),i=1,norb),noa,ydp
      ndegen=ntrans
c
      if(norb.eq.1) go to 99
      ns=nsym
c
610   format(' symetrie par orbitale ',(30i3))
c
c
      read(40) ((isydeg(i,j),i=1,ns),j=1,ns)
      do 1 i=1,ns
c
      do 1 j=1,ns
      its(i,j)=isydeg(i,j)
 612  format(20x,20i3)
1     continue
      do 11 i=1,norb
      if(ntyp.ne.0) go to 11
      itsym(i)=itsyp(i)
      itsyp(i+norb)=itsyp(i)
11    itsym(i+norb)=itsyp(i)
c
c
      if(ntyp.eq.0) ntyp=nsym
      mnorb=norb+norb
c
      do 4 i=1,ntyp
      k=0
      do 4 j=1,mnorb
      if(itsyp(j).ne.i) go to 4
      k=k+1
      ind(j)=k
      indis(k,i)=j
4     continue
      if(ndegen.eq.0) go to 16
      ntrsy=1
      do 15 i=1,ndegen
      read(40 ) j,((isydeg(k,l),k=1,ntyp),l=1,j)
      do 14 l=1,j
      do 14 k=1,mnorb
      inskl=0
      nskl=isydeg(itsyp(k),l)
      ysin=.false.
      if(nskl.ge.0) go to 13
      nskl=-nskl
      ysin=.true.
13    continue
      if(nskl.ne.0) inskl=indis(ind(k),nskl)
      if(ysin) inskl=-inskl
14    isytr(k,ntrsy+l)=inskl
      ntrsy=ntrsy+j
15    continue
16    continue
      igela=igela+1
      if(igelb.lt.norb) igelb=igela+norb
      if(noca.gt.norb) go to 99
      if(.not.yion) go to 109
      mnor1=mnorb+1
      itsym(mnor1)=nsym+1
      nsym1=nsym+1
      its(nsym1,nsym1)=nsym1
      isym=nsym+isym
      do 105 i=1,nsym
c
      its(i,nsym1)=nsym1
105   its(nsym1,i)=nsym1
      iorb(mnor1)=norb+1
      ispin(mnor1)=1
      iwspin(mnor1)=imoins
      isytr(mnor1,1)=0
      if(isz.ne.0) isytr(mnor1,1)=0
      do 107 l=2,ntrsy
      isytr(mnor1,l)=mnor1
107   continue
109   continue
c
      n=norb
      do 150 i=1,norb
      iorb(i)=i
      iorb(i+norb)=i
      ispin(i)=0
      ispin(i+norb)=1
      iwspin(i)=iplus
150   iwspin(i+norb)=imoins
      do 170 i=1,norb
      k=i+norb
      if(k.lt.igelb.and.i.gt.igela) go to 169
      if(k.gt.igelb.and.i.lt.igela) go to 169
      if(k.lt.nocb.and.i.gt.noca) go to 169
      if(k.gt.nocb.and.i.lt.noca) go to 169
      isytr(i,1)=k
      isytr(k,1)=i
      go to 170
169   isytr(i,1)=0
      isytr(k,1)=0
170   continue
      yncper=.false.
      if(ntrsy.lt.2) go to 191
      do 190 n =2,ntrsy
      do 189 k=igela,mnorb
      iskn=isytr(k,n)
      iskn=iabs(iskn)
      if(iskn.lt.igela.or.iskn.gt.mnorb) go to 188
      if(k.gt.noca) go to 185
      if(iskn.gt.noca) go to 188
      go to 189
185   if(k.gt.norb) go to 186
      if(iskn.le.noca.or.iskn.gt.norb) go to 188
      go to 189
186   if(k.gt.nocb) go to 187
      if(iskn.lt.igelb.or.iskn.gt.nocb) go to 188
      go to 189
187   if(iskn.gt.nocb) go to 189
188   isytr(k,n)=0
189   continue
190   continue
191   continue
c
      ncf=0
      ncfpr=0
      ny=0
      yret=.false.
      write(6,*) 'seuil de selection : tau='
      write(6,9999)tau
 9999 format(e12.5,//)
 
 2719 continue
c      write(6,*) 'yret',yret
      call lecnew(ndbl,nts,nps,nexs,nds,yret)

      if(yret) return
      mesp=1
      monot=0
      monop=0
770   continue
c      write(6,*) 'sortie lec',ndbl,(nexs(i),i=1,ndbl)
      if(ndbl.ne.0) then
        do 80 i=1,ndbl
        nec=nexs(i)
        if(nec.eq.0) then
        call transf(0)
        else
c
        ncfo=ncfpr+1
        if(ncfo.gt.ncper) ncfo=ncper+1
        do 20 l=1,nec
        monot=monot+1
        monop=monop+1
        nt(l)=nts(monot)
20      np(l)=nps(monop)
        call transf(nec)
        do 31 l=1,nec
        mt(l)=nt(l)
31      mp(l)=np(l)
        if(ncfpr.lt.ncfo) go to 79
        if(yspipro) then 
        call spipr(maxco,nec,nt,np)
        else 
      if(isz.ne.0) go to 30
      do 25 l=1,nec
      ntl=isytr(nt(l),1)
      npl=isytr(np(l),1)
      if(ntl.eq.0)go to 30
      if(npl.eq.0)go to 30
      nt(l)=ntl
25    np(l)=npl
      call transf(nec)
      end if
30    if(ndegen.eq.0) go to 75
      do 70 ll=2,ntrsy
      ncfa=ncfpr
      isigne=1
      ldeg=0
      do 34 l=1,nec
      nisy=isytr(mt(l),ll)
      if(nisy.lt.0) isigne=-isigne
      n=iabs(nisy)
      if(n.eq.0) go to 70
      nt(l)=n
      if(n.eq.mt(l)) go to 33
      ldeg=ldeg+1
33    nisy=isytr(mp(l),ll)
      if(nisy.lt.0) isigne=-isigne
      n=iabs(nisy)
      if(n.eq.0) go to 70
      if(n.ne.mp(l)) ldeg=ldeg+1
      np(l)=n
34     continue
      if(ldeg.eq.0.or.((2*(ldeg/2)).ne.ldeg))go to 40
      call transf(nec)
      if(ncfpr.le.ncfa) go to 65
      if(yspipro) goto 3400
      goto 3401
 3400 continue
        call spipr(maxco,nec,nt,np)
      goto 3402
 3401 continue
      if(isz.ne.0) go to 65
      do 37 l=1,nec
      ntl=isytr(nt(l),1)
      npl=isytr(np(l),1)
      if(ntl.eq.0) go to 40
      if(npl.eq.0) go to 40
      nt(l)=ntl
37    np(l)=npl
      call transf(nec)
 3402 continue
      if(ncfpr.le.(ncfa+1)) go to 65
40    continue
65    continue
70    continue
75    if(ncfpr.gt.ncper) ncper=ncfo-1
      if(ncfo.ge.ncper) go to 79
      if(ncfpr-ncfo) 79,76,77
76    continue
      go to 79
77    continue
      ncfo=ncfo+1
79    continue
      end if
80    continue
      end if
      if(ncfpr.eq.0) go to 90
      if(ncper.gt.ncfpr) ncper=ncfpr
      if(.true.)go to 85
      yncper=.true.
c
      go to 770
 85   continue
c
c
c
      go to 2719
90    write(6,91)
91    format('  determinant')
      stop
99    write(6,98)
98    format(' erreur dans les informations generales')
9500  write(6,9501)
9501  format(' erreur icinp')
      stop
      end
