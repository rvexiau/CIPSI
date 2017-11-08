      real*8 function hmp(nest,ntst,npst)
      implicit real*8 (a-h,o-x,z),logical*1 (y)
      include 'pshf.prm'
c---------energie moller-plesset du determinant i
      dimension nv1(2),nv2(2)
      equivalence(nv11,nv1(1)),(nv12,nv1(2)),(nv21,nv2(1)),(nv22,nv2(2))
      common/spo/ feff(2*doa),fdiag(2*doa),
     * iorb(2*doa),iwspin(2*doa),ispin(2*doa),
     *itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),isytr(2*doa,nsymz)
      integer*4 iorb,iwspin,itsym,its,isytr
      dimension nest(2),ntst(8,2),npst(8,2)
      common  /cipsi/ norb,mnorb,noca,nocb,
     *igela,igelb,ncf,ncper,nsym,isym,isz,ntrsy,metat,
     *isymat(smax),nelac,nprint,yprt,yion
      common/ic/c(ndetz*metz),e(ndetz),eb(metz),e2(metz),
     *emp(metz),e2en(metz),teste(metz),elevy(metz),e2mp(metz),ac,ei,
     *hdef,dpivot,fmpb(2*doa),dmpb(10),ishif,ietat
      common/truc/indic(doa*(doa+1)/2),jndic(doa*(doa+1)/2),
     1nbo(99),ndeb(500),nad(2000),kt(2000),
     2lndic(doa*(doa+1)/2)
      integer*4 indic,jndic,nad,kt,lndic
      common/jkf/ f(doa*(doa+1)/2),aj(doa*(doa+1)/2),
     *ak(doa*(doa+1)/2),empref(metz),tocom(2*doa)
      num(i)=i*(i-1)/2
      nei=nest(1)
c
      if(nei.eq.0) go to 2300
c
       ac=0.d0
      do 2200 l=1,nei
      n1t=ntst(l,1)
      n1p=npst(l,1)
2200  ac=ac+fmpb(n1p)-fmpb(n1t)
      hmp=ac
      return
2300  hmp=0.d0
      return
      end
