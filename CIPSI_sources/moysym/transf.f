      subroutine transf(nec)
      implicit real*8(a-h,o-x,z),logical*1(y)
      include 'pshf.prm'
      integer*4 trou,part,nexst,nd
      integer*4 iorb,iwspin,itsym,its,isytr
      common/trsf/nt(8),np(8),mt(8),mp(8)
      common  /cipsi/ norb,mnorb,noca,nocb,
     *igela,igelb,ncf,ncper,nsym,isym,isz,ntrsy,metat,
     *isymat(smax),nelac,nprint,yprt,yion
      common /pr/ ncfpr
      common/spo/ feff(2*doa),fdiag(2*doa),
     * iorb(2*doa),iwspin(2*doa),ispin(2*doa),
     *itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),isytr(2*doa,nsymz)
      common/dinde/trou(8*ndimh),part(8*ndimh),nexst(ndimh),
     * nd(ndimh)
      common/e0cip/ebmp(metz),eben(metz),tau

      mnorb=norb+norb
      ncfi=ncfpr+1
      if(nec-1) 10,30,20
10    ncfpr=ncfpr+1
      write(6,100)
100   format(' fondamental')
      call pie(0,nts,nps,ncfpr)
      if(isz.ne.0) ncfpr=ncfpr-1
      isym=1
      ncf=ncfpr
      return
c     ordonnons
c
20    necm=nec-1
      do 25 i=1,necm
      ip=i+1
      do 25 j=ip,nec
      if(nt(j).gt.nt(i)) go to 24
      n=nt(j)
      nt(j)=nt(i)
      nt(i)=n
24    if(np(j).gt.np(i)) go to 25
      n=np(j)
      np(j)=np(i)
      np(i)=n
25    continue
c
30    continue
      is=1
      ls=0
      if(nt(1).lt.igela) go to 80
      if(np(nec).gt.mnorb.and..not.yion)go to 80
      if(nt(nec).gt.nocb) go to 80
      if(np(nec).le. noca) go to 80
      do 35 i=1,nec
      if(nt(i).le.noca) go to 31
      if(nt(i).lt.igelb) go to 80
31    if(np(i).gt.nocb) go to 32
      if(np(i).gt.norb) go to 80
      if(np(i).le.noca) go to 80
32    is=its(is,itsym(nt(i)))
      is=its(is,itsym(np(i)))
      ls=ls+ispin(np(i))
      ls=ls-ispin(nt(i))
35    continue
      if(isym.eq.0) isym=is
      if(.not.yion.and.is.ne.isym) go to 90
      if(ls.ne.isz) go to 96
c    test de non repetition
      do j=1,ncfpr
        jinv=ncfpr-j+1
	if(nexst(jinv).ne.nec)go to 43
	  ndd=nd(jinv)
	  do k=1,nec
	    if(nt(k).ne.trou(ndd+k))go to 43
	    if(np(k).ne.part(ndd+k))go to 43
          end do
c  le determinant est identique a un determinant deja rencontre
	  return
43        continue
      end do
      ncfpr=ncfpr+1
      call pie(nec,nt,np,ncfpr)
      if(nprint.ge.1) write(6,4600) ncfpr,(nt(kk),kk=1,nec),(np(kk),kk
     *=1,nec)
 4600 format(10x,30i4)
      ncf=ncfpr
      return
 46   continue
      if(nprint.ge.1) write(6,4601) ncfpr,(nt(kk),kk=1,nec),(np(kk),kk
     *=1,nec)
 4601 format(i4,' refuse',30i4)
      ncf=ncfpr
      return
70      continue
      ncf=ncfpr
120   return
80    write(6,85)
c
      ncf=ncfpr
      if(yion) return
85    format(' 1 spin orbitale est hors de l espace de base')
      write(6,*)(nt(i),np(i),i=1,nec)
      write(6,*) nec,i,np(nec),np(i),nt(nec),nt(i)
      write(6,*) igela,igelb
      stop
90    write(6,95)
95    format(' determinant de mauvaise symetrie')
      if(nec.ne.0)write(6,*)nec,(nt(j),j=1,nec),(np(j),j=1,nec)
      go to 70
96    write(6,97)
97    format(' mauvais spin')
      write(6,*)(nt(i),np(i),i=1,nec)
      stop
      end
