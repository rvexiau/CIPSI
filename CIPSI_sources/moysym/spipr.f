      subroutine spipr(maxco,nec,nt,np)
      implicit real*8(a-h,o-x,z),logical*1(y)
      include 'pshf.prm'
      integer*4 ist,isp
      common /pr/ ncfpr
      common  /cipsi/ norb,mnorb,noca,nocb,
     *igela,igelb,ncf,ncper,nsym,isym,isz,ntrsy,metat,
     *isymat(smax),nelac,nprint,yprt,yion
      dimension  nt(nexz),np(nexz),isp(doa,nexz+1),ist(doa,nexz+1),
     * mt(nexz),
     * mp(nexz),
     *ycot(nexz+1),ycop(nexz+1),lt(nexz),lp(nexz)
c---------nombre et reperage des couches ouvertes
c
      ncot=nec
      ncop=nec
      if(nprint.ge.2) write(6,*)'det de depart',nec,(nt(i),i=1,nec),
     * (np(i),i=1,nec)
      if(nprint.ge.2) write(6,1000) ncot,ncop
      if(nprint.ge.2) write(6,1000) nt,np,ncfpr,maxco,nec,norb,isz
      do 10 i=1,nec
      ycot(i)=.true.
      ycop(i)=.true.
 10   continue
      do 1 i=1,nec
      if(nt(i).gt.norb) goto 1
      do 2 j=1,nec
      if(nt(j)-norb.ne.nt(i)) goto 2
      ycot(i)=.false.
      ycot(j)=.false.
      ncot=ncot-2
 2    continue
 1    continue
      if(nprint.ge.2) write(6,*) ncot,ycot
      do 3 i=1,nec
      if(np(i).gt.norb) goto 3
      do 4 j=1,nec
      if(np(j)-norb.ne.np(i)) goto 4
      ycop(i)=.false.
      ycop(j)=.false.
      ncop=ncop-2
 4    continue
 3    continue
      if(nprint.ge.2) write(6,*) ncop,ycop
      if(nprint.ge.2) write(6,1000) nt,np,ncot,ncop
 1000 format(30i4)
c---------------stop si trop de couches ouvertes
      nco=ncot+ncop
      if(nco.le.maxco) goto 5
      if(nprint.ge.0) write(6,1001) (nt(k),k=1,nec),(np(k),k=1,nec)
 1001 format(' trop de couches ouvertes sur le determinant',12i4)
      stop
  5   continue
      if(nco.le.8) goto 6
      if(nprint.ge.0) write(6,1001) (nt(k),k=1,nec),(np(k),k=1,nec)
      stop
 6    continue
c-------------on fabrique un determinant de depart tous spins ^
c-------------sur les couches ouvertes
      do 11 i=1,nec
      mt(i)=nt(i)
      mp(i)=np(i)
      lt(i)=nt(i)
      lp(i)=np(i)
      if(.not.ycot(i)) goto 12
      if(mt(i).gt.norb) mt(i)=mt(i)-norb
 12   continue
      if(.not.ycop(i)) goto 11
      if(mp(i).gt.norb) mp(i)=mp(i)-norb
 11   continue
      if(nprint.ge.1) write(6,*)'nec,mt,mp',nec,(mt(k),k=1,nec),
     * (mp(k),k=1,nec)
 1002 format(30x,20i4)
      if(nprint.ge.2) write(6,1000) mt,mp,lt,lp
c------------------degres de liberte sur le choix des + et des -
c------------------des trous (id. pour les particules)
      minco=min0(ncot,ncop)
      maxnpt=(ncot+minco)/2+1+isz
      minnpt=(ncot-minco+1)/2+1+isz
      maxnpp=(ncop+minco)/2+1
      if(maxnpt-1.gt.ncot) maxnpt=ncot+1
      if(isz.ne.0) minnpt=minnpt-(ncop-minco)/2
      if (minnpt.le.0) minnpt=1
      if(minnpt.gt.maxnpt) minnpt=maxnpt
      if(nprint.ge.2) write(6,*) minco,minnpt,maxnpt,maxnpp
c
c------------------------generation des equiv de spin
c
c------------boucle sur le nb de trous possible
      do 13 it1=minnpt,maxnpt
      it=it1-1
      if(nprint.ge.2) write(6,1000) it
      if(maxnpt.eq.1) goto 33
      if(it.eq.0) goto 33
      if(it.eq.ncot) goto 33
      call combin(ist,ncot,it,ntt,doa,nexz+1)
      if(nprint.ge.2) write(6,1000) ncot,it,ntt
      if(nprint.ge.2) write(6,*)((ist(i5,j5),j5=1,it),i5=1,ntt)
      goto 34
 33   ntt=1
 34   continue
      do 14 jt=1,ntt
      if(maxnpt.eq.1) goto 30
      if(it.eq.0) goto 20
      if(it.eq.ncot) goto 21
      jj=1
      jco=0
      do 15 k=1,nec
      nt(k)=mt(k)+norb
      if(ycot(k)) goto 152
      nt(k)=nt(k)-norb
      goto 15
 152  continue
      if(jj.gt.it) goto 15
      jco=jco+1
      if(ist(jt,jj).ne.jco) goto 15
      nt(k)=nt(k)-norb
      jj=jj+1
 15   continue
      if(nprint.ge.2) write(6,1000) nt
      goto 23
 20   do 22 i=1,nec
      nt(i)=mt(i)
 22   if(ycot(i)) nt(i)=mt(i)+norb
      goto 23
 21   do 24 i=1,nec
 24   nt(i)=mt(i)
      goto 23
 30   do 31 i=1,nec
 31   nt(i)=lt(i)
 23   continue
      if(nprint.ge.2) write(6,*) nt
c-----------------on a genere les trous d'un determinant
c------------boucle sur les particules (interne)
      ip=it-isz+(ncop-ncot)/2
      if(nprint.ge.2) write(6,*) ip
      if(maxnpp.eq.1) goto 35
      if(ip.eq.0) goto 35
      if(ip.eq.ncop) goto 35
      call combin(isp,ncop,ip,ntp,doa,nexz+1)
      goto 36
 35   ntp=1
 36   continue
      do 17 jp=1,ntp
      if(maxnpp.eq.1) goto 32
      if(ip.eq.0) goto 25
      if(ip.eq.ncop) goto 26
      jj=1
      jco=0
      do 18 k=1,nec
      np(k)=mp(k)+norb
      if(ycop(k)) goto 182
      np(k)=np(k)-norb
      goto 18
 182  continue
      if(jj.gt.ip) goto 18
      jco=jco+1
      if(isp(jp,jj).ne.jco) goto 18
      np(k)=np(k)-norb
      jj=jj+1
 18   continue
      goto 28
 25   do 27 i=1,nec
      np(i)=mp(i)
 27   if(ycop(i)) np(i)=mp(i)+norb
      goto 28
 26   do 29 i=1,nec
 29   np(i)=mp(i)
      goto 28
 32   do 37 i=1,nec
 37   np(i)=lp(i)
 28   continue
      if(nprint.ge.2) write(6,*) 'det genere' , (nt(k),k=1,nec),
     * (np(k),k=1,nec),ncfpr
      if(yion.and.np(nec).ne.norb*2+1) goto 17
      call transf(nec)
 17   continue
 14   continue
 13   continue
      return
      end
