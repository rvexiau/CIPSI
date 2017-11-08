      subroutine morue
      implicit real*8(a-h,o-x,z),logical*1(y)
      include 'pshf.prm'
      integer*4 trou,part,nexst,nd
      integer*4 ibid
      character*20 title
      character*40 typener
      common/vect/v(ndimh),stefmp(ndimh),stefen(ndimh)
      common/dinde/trou(8*ndimh),part(8*ndimh),nexst(ndimh),
     * nd(ndimh)
      common/e0cip/ebmp(metz),eben(metz),tau
      common/ic/c(ndetz*metz),e(ndetz),eb(metz),e2(metz),
     *emp(metz),e2en(metz),teste(metz),elevy(metz),e2mp(metz),ac,ei,
     *hdef,dpivot,fmpb(2*doa),dmpb(10),ishif,ietat
      common  /cipsi/ norb,mnorb,noca,nocb,
     *igela,igelb,ncf,ncper,nsym,isym,isz,ntrsy,metat,
     *yprt,yion,ystkmp
      common/nom33/title
      dimension  cunt(metz)
      data zero/0.d0/
      dimension numero(metz),hefmp(nhefz),
     *hefen(nhefz),hefenb(nhefz),
     *hef(nhefz)
      integer itab(lsize),jtab(lsize)
      dimension vbuf(lsize)
      equivalence (v(1),vbuf(1))
c
c
      read(4) norb,noca,nocb,metat,ndims,isz,yion,ii,jj,
     * (ibid,ibid,i=1,ndims),
     & (ibid,ibid,i=1,ii),(c(i),i=1,jj),(e(i),i=1,metat),
     &  (ebmp(i),i=1,metat),(eben(i),i=1,metat),(numero(i),i=1,metat),
     &   (xdum,i=1,ndims),(fmpb(i),i=1,norb)
      mmm=metat*(metat+1)/2
      read(4)
      read(4) (hef(i),i=1,mmm)
      
c      write(6,9342)
 9342 format(//,' elements diagonaux m.p')
c       write(6,9343) (stefmp(i),i=1,ny)
 9343 format(12f10.6)
      write(6,9220)
9220  format(//)
      write(6,9344)
 9344 format(/' caracteristique des etats a l''ordre zero')
      write(6,9230)
9230  format(1x,40('*'),/)
      write(6,9233)
9233  format('  etat',4x,'partition m.p.b',10x,'partition e.n.vp',10x,
     *'partition e.n.b',/)
      do 9345 i=1,metat
 9345 write(6,9346) i,ebmp(i),e(i),eben(i)
 9346 format(i4,f18.8,2f25.8)
      ij=0
      do 25 i=1,metat
      cunt(i)=0.d0
      do 25 j=1,i
      ij=ij+1
       hefenb(ij)=0.d0
       hefmp(ij)=0.d0
25     hefen(ij)=0.d0
      last=lsize
      itab(lsize)=0
      nbuf=lsize
      if(ndims.eq.ncf) go to 8300
    2 nbuf=nbuf+1
      if(nbuf.le.last) go to 5
      if(itab(lsize).eq.-1) go to 4
      read(1) itab,jtab,vbuf
      if(itab(lsize).eq.-1) last=jtab(lsize)
      nbuf=1
 5    k=itab(nbuf)
      i=jtab(nbuf)
c    k numero du determinant exterieur
c    i numero du determinant dans s
      if(k.eq.1) go to 2
c     k=1  element diagonal a sauter
      if(k.eq.i.and.k.gt.ndims) go to 7
c     k=i   dernier terme pour k donne on peut accumuler les cont.
c     dans hefmp,etc
      if(k.le.ndims.or.i.gt.ndims) go to 2
c     pour que l'element de matrice contribue il faut que
c     k exterieur a s et i dans s
      hij=vbuf(nbuf)
      mi=i-ndims
      do 102 m=1,metat
      mi=mi+ndims
      cunt(m)=cunt(m)+c(mi)*hij
  102 continue
      go to 2
7     stenk=stefen(k)
      stmpk=stefmp(k)
      enk2=stenk+stenk
      empk2=stmpk+stmpk
      ii=0
      do 103 m=1,metat
c     write(6,*) k,stmpk,stenk,cunt(1),eben(1),emp(1),e(1)
      do 103 n=1,m
      ii=ii+1
      denb=eben(m)+eben(n)-enk2
      dmp=ebmp(m)+ebmp(n)-empk2
      den=e(m)+e(n)-enk2
      hmn=cunt(m)*cunt(n)
      hmn=hmn+hmn
      hefmp(ii)=hefmp(ii)+hmn/dmp
      hefen(ii)=hefen(ii)+hmn/den
      hefenb(ii)=hefenb(ii)+hmn/denb
103   continue
      do 105 m=1,metat
  105 cunt(m)=0.d0
      go to 2
4     continue
      write(6,9000)
9000  format(//,1x,'contribution 2n ordre des determinants moyens',/)
      write(6,9234)
9234  format(1x,45('*'))
      write(6,9001)
9001  format(1x,'perturbation mp',/)
      do 120 m=1,metat
      i1=m*(m-1)/2+1
      i2=m*(m+1)/2
120   write(6,9004) (hefmp(ii),ii=i1,i2)
      write(6,9220)
      write(6,9003)
9003  format(1x,'perturbation envp',/)
      do 140 m=1,metat
      i1=m*(m-1)/2+1
      i2=m*(m+1)/2
140   write(6,9004) (hefen(ii),ii=i1,i2)
      write(6,9220)
      write(6,9002)
9002  format(1x,'perturbation enb',/)
      do 130 m=1,metat
      i1=m*(m-1)/2+1
      i2=m*(m+1)/2
130   write(6,9004) (hefenb(ii),ii=i1,i2)
9004  format(11f11.6)
      write(6,9220)
      do m=1,metat
	cunt(m)=hefmp(m*(m+1)/2)
      end do 
      typener='moyens mbp '//title
      call stkener(typener,cunt,metat)
      do m=1,metat
	cunt(m)=hefen(m*(m+1)/2)
      end do 
      typener='moyens envp'//title
      call stkener(typener,cunt,metat)
      do m=1,metat
	cunt(m)=hefenb(m*(m+1)/2)
      end do 
      typener='moyens enb '//title
      call stkener(typener,cunt,metat)
      write(6,9005)
      write(6,9005)
9005  format(/,2x,128('*'))
 8300 rewind 4
      read(4)
      read(4)
      read(4)
      write(4) (hefmp(i),i=1,mmm),
     *         (hefen(i),i=1,mmm),
     *         (hefenb(i),i=1,mmm)
     *        ,ny
      rewind 1
      return
      end
