      real*8 function hntd(nest,ntst,npst)
      implicit real*8 (a-h,o-x,z),logical*1 (y)
      include 'pshf.prm'
c        but
c
c                  ssp de recherche des elements d interaction
c                  entre deux determinants,numerote ii  et  jj
c
      dimension ita(20),itb(20),nv1(2),nv2(2),natu1(2),natu2(2)
      equivalence(nv11,nv1(1)),(nv12,nv1(2)),(nv21,nv2(1)),(nv22,nv2(2))
      common/spo/ feff(2*doa),fdiag(2*doa),
     * iorb(2*doa),iwspin(2*doa),ispin(2*doa),
     *itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),isytr(2*doa,nsymz)
      integer*4 iorb,iwspin,itsym,its,isytr
      dimension nest(2),ntst(8,2),npst(8,2)
      common  /cipsi/ norb,mnorb,noca,nocb,
     *igela,igelb,ncf,ncper,nsym,isym,isz,ntrsy,metat,
     *yprt,yion
      common/ic/c(ndetz*metz),e(ndetz),eb(metz),e2(metz),
     *emp(metz),e2en(metz),teste(metz),elevy(metz),e2mp(metz),ac,ei,
     *hdef,dpivot,fmpb(2*doa),dmpb(10),ishif,ietat
      common/truc/indic(doa*(doa+1)/2),jndic(doa*(doa+1)/2),
     1nbo(99),ndeb(500),nad(2000),kt(2000),
     2lndic(doa*(doa+1)/2)
      integer*4 indic,jndic,nad,kt,lndic
      common/jkf/ f(doa*(doa+1)/2),aj(doa*(doa+1)/2),
     *ak(doa*(doa+1)/2),empref(metz),tocom(2*doa)
      logical iden
      num(i)=i*(i-1)/2
c      write(6,*) 'hntd', nest,ntst,npst
      i=1
      j=2
      iden=.false.
c
c                  le determinant - 1 - est choisi comme le
c                  determinant ayant le plus d o.m. excitees
c
      if (nest(i)-nest(j)) 5,5,7
 7    kr = i
      kl= j
      go to 9
 5    kr= j
      kl= i
 9    ne1=nest(kr)
c                  test sur le fondamental
      if(ne1.eq.0) goto 44
      ne2=nest(kl)
c                  nbdif : est le nombre de spin-orbitales
c                          differentes
      nbdif=ne1-ne2
c
      if(nbdif.gt.2) goto 24
c                  construction pour le determinant -1- et -2-
c                  de it  et ns,qui contiennent les numeros
c                  des orbitales et leurs spins
c
c
c
      do 15 j=1,ne1
      ita(j)=ntst(j,kr)
      ita(j+ne1)=npst(j,kr)
15    continue
c
c                  ncr :  est le nombre de croisements
c                         dans le diagramme d interaction
      ysig=.false.
c                  si ne2=0 le determinant -2- est le fondamental
      if (ne2.eq.0) goto 125
121   continue
      do 25 j= 1,ne2
      itb(j)=ntst(j,kl)
 25   itb(j+ne2)=npst(j,kl)
    4 do 8 k=1,2
      k1=(k-1)*ne1
      k2=(k-1)*ne2
      do 10i=1,ne2
c
      nit1=itb(k2+i)
      j1=1+k1
      j2=ne1+k1
      do 12 j=j1,j2
      njt2=ita(j)
      if(nit1.ne.njt2) go to 12
      if(i.eq.(j-k1)) go to 10
      ita(j)=ita(i+k1)
      ita(i+k1)=njt2
      ysig=.not.ysig
      go to 10
   12 continue
      nbdif=nbdif+1
      if(nbdif.gt.2) go to 24
      nv1(nbdif)=nit1
c
      natu1(nbdif)=k
      natu2(nbdif)=i
   10 continue
    8 continue
125   if(nbdif.le.0) go to 44
c
c                   nar : excitation suplementaire du determinant -1-
c                        par rapport au determinant -2-
   26 nar=ne1-ne2
      if(nar.le.0) go to 28
   30 do 32 i=1,nar
      net=ne1+1-i
      nv1(i)=ita(net)
      nv2(i)=ita(net+ne1)
   32 continue
   28 nbar=nbdif-nar
      if(nbar.le.0) go to 34
   36 nar1=nar+1
      do 38 i=nar1,nbdif
      k=natu1(i)
      ni=natu2(i)
      nu=nv1(i)
      nv2(i)=ita(ni+ne1*(k-1))
c
      if(k.eq.2) go to 38
40    nv1(i)=nv2(i)
      nv2(i)=nu
      ysig=.not.ysig
   38 continue
34    continue
c                  calcul de l element de matrice
 2    ha=0.
      n1=iorb(nv11)
      n2=iorb(nv21)
      ns1=ispin(nv11)
      ns2=ispin(nv21)
      ncr=1
      if(ysig) ncr=-1
      if(nbdif.eq.1) go to 58
c                  les -2- determinants differe par -2- spin-orbitales
60    n3=iorb(nv12)
      n4=iorb(nv22)
      ns3=ispin(nv12)
      ns4=ispin(nv22)
      if(ns1.eq.ns2)ha=ai(n1,n2,n3,n4)
      if(ns1.eq.ns4)ha=ha-ai(n1,n4,n3,n2)
      go to 75
c
c                  les -2- determinants differe par -1- spin-orbitale
58    if(n2.ge.n1) go to 57
      n12=num(n1)+n2
      go to 56
57    n12=num(n2)+n1
56    ha=f(n12)
      do 72 i=1,ne1
      np=iorb(ita(i))
      nsu=ispin(ita(i))
      ha=ha-ai(n1,n2,np,np)
      if(ns1.eq.nsu) ha=ha+ai(n1,np,n2,np)
      np=iorb(ita(i+ne1))
      if(np.gt.norb) go to 72
      nsu=ispin(ita(i+ne1))
      ha=ha+ai(n1,n2,np,np)
      if(ns1.eq.nsu) ha=ha-ai(n1,np,n2,np)
   72 continue
75    hntd=ha*ncr
      if(.not.yprt) return
      write(6,1007)hntd
 1007 format(10x,f15.10)
      return
 24   hntd=0.0
      if(.not.yprt) return
      write(6,1005)
 1005 format(1h0,'pas d''interaction')
      return
   44 iden=.true.
      hntd=0.0
      if(.not.yprt) return
      write(6,1006)
 1006 format (1h0,'identite')
      return
      end
