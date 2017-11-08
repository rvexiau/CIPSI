      program moy
      implicit real*8(a-h,o-x,z),logical*1(y)
      include 'pshf.prm'
      character*20 title
      character*1 jspin
!      character(len=20) fname
      integer*4 trou,part,nexst,nd
      integer*4 jdet
      integer*4 iorb,iwspin,itsym,its,isytr,kkk,lll
      integer   itab(lsize),jtab(lsize)
      common/spo/ feff(2*doa),fdiag(2*doa),
     * iorb(2*doa),iwspin(2*doa),ispin(2*doa),
     *itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),isytr(2*doa,nsymz)
      common  /cipsi/ norb,mnorb,noca,nocb,
     *igela,igelb,ncf,ncper,nsym,isym,isz,ntrsy,metat,
     *yprt,yion,ystkmp
      common/ic/c(ndetz*metz),e(ndetz),eb(metz),e2(metz),
     * emp(metz),e2en(metz),teste(metz),elevy(metz),e2mp(metz),ac,ei,
     *hdef,dpivot,fmpb(2*doa),dmpb(10),ishif,ietat
      common/dinde/trou(8*ndimh),part(8*ndimh),nexst(ndimh),
     * nd(ndimh)
      common/mod90/noa,ydp
      common/e0cip/ebmp(metz),eben(metz),tau
      common/nom33/title
      common/vect/v(ndimh),stefmp(ndimh),stefen(ndimh)
      dimension ntst(8,2),npst(8,2),nest(2)
      dimension vbuf(lsize)
      dimension jspin(12),jdet(12)
      namelist/moyen/yprt,tau,yecr,ystkmp,ymat,iord,title,yrcmp,
     * y_secd,yspipro,ysymat
c
      call openf
      trt=0
 0901 tau=0.d0
      ymat=.false.
      yprt=.false.
      yecr=.false.
      ystkmp=.false.
      yrcmp=.false.
      y_secd=.false.
      yspipro=.false.
      title='                    '
      numdet=0
      nhij=0
      nhntd=0
      rewind 4
      write(6,*) 'bada bada'
      read(4,end=9000,err=9900)norb,noca,nocb,metat,ncf,isz,yion,ii,jj,
     y (kkk,lll,i=1,ncf),
     & (kkk,lll,i=1,ii),(xdum,i=1,jj),(dumx,dumy,dumz,m=1,metat),
     & (jdum,i=1,metat),(xdum,i=1,ncf),(fmpb(i),i=1,norb)
      write(6,*) 'bada bada'
      ndet_s=ncf
      rewind 4
      mnorb=norb+norb
      do 2567 i=1,norb
 2567 fmpb(i+norb)=fmpb(i)
      fmpb(norb+norb+1)=0.d0
      read(5,moyen)
      write(6,7854) tau,ndimh,yspipro
 7854 format(
     * ' seuil de selection ',f10.5, /,
     * ' nombre maximal de determinants ',i10,/,
     * ' generation de toutes les combinaisons de spin(yspipro) ',l4)
      call deter(yspipro)
      write(6,160) 
  160 format(' determinants engendres')
      write(6,*) ncf
      if(ncf.gt.ndimh) stop 99
111   format(1x,30i4)
c   si ypsipro=t, on cree une file 60-bis (62) contentant des fonctions 
c propres de spin
c
      if(yspipro) call stk62(ncf)
!      write(6,*)'ATTTEEEENTION: ncf=',ncf
!      inquire(61,name=fname)
!      write(6,*)'moy writing to ',fname

c     RV2015 : informations sur les determinants pour la file 64
      ii=nd(ncf)+nexst(ncf)+1
      write(64)ncf,ii,mnorb,yion,(nexst(i),nd(i),i=1,ncf),
     *(trou(i),part(i),i=1,ii),
     *(iorb(i),i=1,mnorb+1),(ispin(i),i=1,mnorb+1) 
      do i=1,ncf
      ne1=nexst(i)
      write(61)ne1
      !write(6,*)'ne1=',ne1
      if(ne1.ne.0)then
	write(61)(trou(nd(i)+j),j=1,ne1),(part(nd(i)+j),j=1,ne1)
!	write(6,*)(trou(nd(i)+j),j=1,ne1),(part(nd(i)+j),j=1,ne1)
      endif
      if(yecr) then 
      do  k=1,ne1
        jspin(k)='+'
        jspin(k+ne1)='+'
        if(trou(nd(i)+k).gt.norb) jspin(k)='-'
        if(part(nd(i)+k).gt.norb) jspin(k+ne1)='-'
        jdet(k)=trou(nd(i)+k)
        jdet(k+ne1)=part(nd(i)+k)
        if(jdet(k).gt.norb) jdet(k)=jdet(k)-norb
        if(jdet(k+ne1).gt.norb) jdet(k+ne1)=jdet(k+ne1)-norb
      end do
      write(6,9998)i, (jdet(k),jspin(k),jdet(k+ne1),jspin(k+ne1),
     *k=1,ne1)
 9998 format(1x,i5,3x,20(i3,a1,1x))
      end if
      end do
c si yrcmp=true, on ne calcule pas la matrice h donc fin du programme
      if(yrcmp)stop
c
c
      call ijkf
      call reijkl(norb,nsym,ydp)
      my1=0
      nbuf=0
      do 3000 i=1,mnorb
 3000 yoc(i)=.false.
      yoc(mnorb+1)=.false.
c     boucle sur les determinants
      do 1 i=1,ncf
      ne1=nexst(i)
      if(ne1.eq.0) go to 100
      ndd=nd(i)
      do k=1,ne1
      ntk=trou(ndd+k)
      ntst(k,1)=ntk
      yoc(ntk)=.true.
      npk=part(ndd+k)
      npst(k,1)=npk
      yoc(npk)=.true.
      end do
 100  nest(1)=ne1
        my2=0
        im1=i-1
        if(im1.eq.0) go to 7
	if(y_secd.and.im1.gt.ndet_s) im1=ndet_s
        do 6 j=1,im1
          ne2=nexst(j)
          ndif=ne1-ne2
          if(ndif.gt.2) go to 6
          if(ne2.eq.0) go to 101
          ndd=nd(j)
          do k=1,ne2
            ntk=trou(ndd+k)
            ntst(k,2)=ntk
            if(.not.yoc(ntk))ndif=ndif+1
            npk=part(ndd+k)
            npst(k,2)=npk
            if(.not.yoc(npk))ndif=ndif+1
          end do
          if(ndif.gt.2) go to 6
 101      nest(2)=ne2
          hij=hntd(nest,ntst,npst)
          nhntd=nhntd+1
          if(dabs(hij).lt.1.d-10) go to 6
          nhij=nhij+1
          nbuf=nbuf+1
          itab(nbuf)=i
          jtab(nbuf)=j
          vbuf(nbuf)=hij
          if(nbuf.lt.lsize) go to 6
          write(1) itab,jtab,vbuf
 9547     format(7(2i4,d10.2))
          if(ymat) write(6,9547) (itab(k),jtab(k),vbuf(k),k=1,lsize)
          nbuf=0
 6      continue
    7   stfn=hdig(nest,ntst,npst)
        nhntd=nhntd+1
        nhij=nhij+1
        stefen(i)=stfn
        stefmp(i)=hmp(nest,ntst,npst)
        nbuf=nbuf+1
        itab(nbuf)=i
        jtab(nbuf)=i
        vbuf(nbuf)=stfn*0.5d0
        if(nbuf.lt.lsize) go to 8
        write(1) itab,jtab,vbuf
        if(ymat) write(6,9547) (itab(k),jtab(k),vbuf(k),k=1,lsize)
        nbuf=0
 8      if(ne1.eq.0) go to 1
        do 3010 k=1,ne1
          ntk=ntst(k,1)
          yoc(ntk)=.false.
          npk=npst(k,1)
          yoc(npk)=.false.
 3010   continue
 1      continue
        itab(lsize)=-1
        jtab(lsize)=nbuf
        write(1) itab,jtab,vbuf
        rewind 1
        if(ymat) write(6,9547) (itab(k),jtab(k),vbuf(k),k=1,lsize)
        if(nbuf.eq.0) write(6,845)
  845   format(//' **** attention le dernier buffer sur 1 ne contient'
     y 'aucun element')
        write(6,846) nhntd,nhij
  846   format(//' *** fin de la construction de h ***',/,
     y ' appels a hntd ',i10,' elements stockes sur 1 ',i10)
c      write(6,3456)
c      write(6,3457) (stefen(i),i=1,ncf)
 3456 format(//' elements diagonaux')
 3457 format(12f10.6)
      if(ystkmp) then 
        write(63) (stefmp(i),i=1,ncf)
        write(6,3456)
        do i=1,ncf
         write(6,*) i,stefen(i),stefmp(i)
        enddo
        write(6,2378)
 2378 format(/, ' elements diagonaux m.p. sur file 63')
      endif
      
       call morue
!       if (trt.eq.0) then
!         rewind 1
!         rewind 4
!         rewind 5
!         rewind 40
!         rewind 60
!         rewind 61
!         rewind 62
!         rewind 33
!         trt=1
!         goto 0901
!       endif   
 9000 write(6,9010)
 9010 format(//, '********* fin de fichier sur 4   ******'
     y 'programme stoppe')
      stop
 9900 write(6,9910)
 9910 format(//, '********* erreur fichier sur 4   ******'
     y 'programme stoppe')     
      end
