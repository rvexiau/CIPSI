      program moy
      implicit real*8(a-h,o-x,z),logical*1(y)
      include 'pshf.prm'
      character*20 title
      character*1 jspin
!      character(len=20) fname
      integer*4 trou,part,nexst,nd
      integer*4 jdet
      integer*4 iorb,iwspin,itsym,its,isytr,kkk,lll
      common/spo/ feff(2*doa),fdiag(2*doa),
     * iorb(2*doa),iwspin(2*doa),ispin(2*doa),
     *itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),isytr(2*doa,nsymz)
      common  /cipsi/ norb,mnorb,noca,nocb,
     *igela,igelb,ncf,ncper,nsym,isym,isz,ntrsy,metat,
     *isymat(smax),nelac,nprint,yprt,yion,ystkmp
      common/ic/c(ndetz*metz),e(ndetz),eb(metz),e2(metz),
     * emp(metz),e2en(metz),teste(metz),elevy(metz),e2mp(metz),ac,ei,
     *hdef,dpivot,fmpb(2*doa),dmpb(10),ishif,ietat
      common/dinde/trou(8*ndimh),part(8*ndimh),nexst(ndimh),
     * nd(ndimh)
      common/mod90/noa,ydp
      common/e0cip/ebmp(metz),eben(metz),tau
      common/nom33/title
      dimension ntst(8,2),npst(8,2),nest(2)
      dimension jspin(12),jdet(12)
      namelist/moyen/yprt,tau,yecr,ystkmp,ymat,iord,title,yrcmp,
     * y_secd,yspipro,isymat,nelac,nprint
c
    
      call openf
      trt=0
 0901 tau=0.d0
      isymat=1
      ymat=.false.
      yprt=.false.
      yecr=.false.
      ystkmp=.false.
      yrcmp=.false.
      y_secd=.false.
      yspipro=.false.
      nprint=0
      title='                    '
      numdet=0
      nelac=0
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
c      ecriture de la matrice h

       call morue


 9000 write(6,9010)
 9010 format(//, '********* fin de fichier sur 4   ******'
     y 'programme stoppe')
      stop
 9900 write(6,9910)
 9910 format(//, '********* erreur fichier sur 4   ******'
     y 'programme stoppe')     
      end
