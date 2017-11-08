      real*8 function hdig(nest,ntst,npst)
      implicit real*8 (a-h,o-x,z),logical*1 (y)
      include 'pshf.prm'
c------calcul de l'element diagonal de h
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
      common/jkf/ f(doa*(doa+1)/2),aj(doa*(doa+1)/2),
     *ak(doa*(doa+1)/2),empref(metz),tocom(2*doa)
      num(i)=i*(i-1)/2
      nei=nest(1)
      if(nei.eq.0) go to 1300
      ab=0.
      do 1200 l=1,nei
      n1t=ntst(l,1)
      n1p=npst(l,1)
      ab=ab+fdiag(n1p)-fdiag(n1t)
      ns1t=ispin(n1t)
      ns1p=ispin(n1p)
      n1t=iorb(n1t)
      n1p=iorb(n1p)
c
c
      do 1100 j=1,l
      n2t=ntst(j,1)
      n2p=npst(j,1)
      ns2t=ispin(n2t)
      ns2p=ispin(n2p)
      n2t=iorb(n2t)
      n2p=iorb(n2p)
      if(n2t.gt.n1t) go to 1010
      n12t=num(n1t) +n2t
      go to 1020
1010  n12t=num(n2t)+n1t
1020  if(n2p.gt.n1p) go to 1030
      n12p=num(n1p)+n2p
      go to 1040
1030  n12p=num( n2p)+n1p
1040  n1tp=num(n2p)+n1t
      n2tp=num(n1p)+n2t
      ca=aj(n12t)+aj(n12p)-aj(n1tp)-aj(n2tp)
      if(ns1t.eq.ns2t) ca=ca-ak(n12t)
      if(ns1p.eq.ns2p) ca=ca-ak(n12p)
      if(ns1t.eq.ns2p) ca=ca+ak(n1tp)
      if(ns1p.eq.ns2t) ca=ca+ak(n2tp)
      if(j.eq.l) ca=0.5*ca
1100  ab=ab+ca
1200  continue
      hdig=ab
c
      return
1300  hdig=0.d0
c
      return
      end
