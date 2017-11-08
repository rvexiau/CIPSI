      subroutine itijkl(bijkls,bijkld,ydp)
c       ipair,nad,ktypsy static size increased      
      implicit real*8 (a-h,o-x,z),logical*1 (y)
      include 'pshf.prm'    
c     ydp= t  calcul double precision
c     ydp= f  calcul simple precision 
      real*8,dimension(:),allocatable :: xx_pqrs
      integer,dimension(:),allocatable :: ix_pqrs     
      integer*4 wijkl
      parameter (idinir=72)
      real*4 bijkls,bufs1,bufs2
      common norb,nsym,its(nsymz,nsymz),itsym(doa),itsyv(2*doa),
     1 nvec,isymp(nsymz),mijkl,nrec,mblk
      common isydeg(nsymz,nsymz),ntrsy,istop,iarwr(2*doa),
     1 nocsy(nsymz),ymono,yprt,yprtv,yatom,ywvect,ynocsy
      common /om/ c(doa*doa*2), nbp(nsymz),nbo(nsymz),ind(nsymz+1),
     y ivecsy(2*doa*nsymz)
      common /hondo/ nt,ishell(doa),kloc(doa),ktyp(doa),
     1 numc(doa),newsh(doa,nsymz),iptr(3,nsymz),idtr(6,nsymz),
     2 iftr(10,nsymz),igtr(15,nsymz),yhondo,yprth
       common/tampon/bufs1(lsize),bufs2(lsize),bufd1(lsize),
     1 bufd2(lsize),monowr

       dimension ipair(2000),nad(10000)
     *  ,ivnir(idinir,9),yvnir( idinir,3)
      dimension bijkls(*),bijkld(*)
      integer*2 nad,ktypsy(10000),ltabs(10,10)
       dataltabs/ 9*1,10,
     &                1,2,3,4,5,6,7,8,9,1,
     &                1,3,2,5,4,7,6,9,8,1,
     &                1,4,5,2,3,8,9,6,7,1,
     &                1,5,4,3,2,9,8,7,6,1,
     &                1,6,7,8,9,2,3,4,5,1,
     &                1,7,6,9,8,3,2,5,4,1,
     &                1,8,9,6,7,4,5,2,3,1,
     &                1,9,8,7,6,5,4,3,2,1,10,9*1/
c
      numnu(i)=(i*(i-1))/2
      mrec=nrec
      mblk1=mblk         

      inijkl=lsize
      nsyt=nsym*(nsym+1)/2
      nsyt=nsyt*(nsyt+1)/2
      nsyt=min0(2000,nsyt)
      do 30 i=1,nsyt
      ktypsy(i)=1
30    nad(i)=0
      npaire=0
      npair=npaire+npaire
      nsymp=isymp(nsym)
            ii=0
      nrec=0
      mijkl=0
      if(ymono) go to 690
        ncpqkl=0         
c     RV 01/2016 : store pqrs file in memory
      imin=1
      imax=lsize  
      rewind 8
      allocate(xx_pqrs(lsize*(mrec+1)),ix_pqrs(lsize*(mrec+1)))
      do i=1,mrec     
      read(8)xx_pqrs(imin:imax),ix_pqrs(imin:imax)     
      imin=imin+lsize
      imax=imin+lsize-1
      enddo
      read(8,end=100)xx_pqrs(imin),ix_pqrs(imin)
100   continue

      do 600 msi=1,nsymp
      do 600 msj=1,nsymp
c
      ypq=msi.eq.msj
      do 600 msk=1,nsymp
      mslm=msk
      if(msk.eq.msi) mslm=msj
      do 600 msl=1,nsymp
      yr=.true.
      yrs=msk.eq.msl
      ypqrs=msi.eq.msk.and.msj.eq.msl
          ncpqkl=0
      nrec=nrec+1
      do 500 ii=1,nsym
c
      do 500 jj=1,ii
c
      mono1=nbp(ii)*nbp(jj)
      if(ypq) mono1=(mono1+nbp(jj))/2
      yij=ii.eq.jj
      nij=ii*(ii-1)/2+jj
      ijs=its(ii,jj)
      mij=nbo(ii)*nbo(jj)
      if(yij) mij=(mij+nbo(ii))/2
      do 500 kk=1,ii
c
      llm=kk
      if(kk.eq.ii) llm=jj
      do 500 ll=1,llm  
      mono2=nbp(kk)*nbp(ll)
      if(yrs) mono2=(mono2+nbp(ll))/2
c      write(6,*) 'ii,jj,kk,ll,isymp,mono1,mono2',
c     * ii,jj,kk,ll,isymp(ii),isymp(jj),isymp(kk),isymp(ll),mono1,mono2
      if(isymp(ii).ne.msi.or.isymp(jj).ne.msj.or.
     & isymp(kk).ne.msk.or.isymp(ll).ne.msl) go to 500
      if(mono1.le.0.or.mono2.le.0) go to 500
      ykl=kk.eq.ll
      nijkl=mono1*mono2
      if(ypqrs) nijkl=(mono1+nijkl)/2
      nijkl1=nijkl+1
      yijkl=ii.eq.kk.and.jj.eq.ll
c
      kls=its(kk,ll)
      if(its(kls,ijs).ne.1) go to 500
      nkl=kk*(kk-1)/2+ll
      mkl=nbo(kk)*nbo(ll)
      if(ykl) mkl=(mkl+nbo(kk))/2
      ijkls=nij*(nij-1)/2+nkl
      wijkl=mkl*mij
      if(wijkl.le.0) go to 500
      if(yijkl) wijkl=(wijkl+mij)/2
       if(ntrsy.eq.0)  go to 205
      do 200 n=1,ntrsy
      nklst=nkl
      isig=1
         iit=isydeg(ii,n)
      if(iit) 170,200,171
170     iit=-iit
       isig=-1
171      jjt=isydeg(jj,n)
      if(jjt) 172,200,173
172     jjt=-jjt
      isig=-isig
173      kkt=isydeg(kk,n)
      if(kkt) 174,200,175
174     kkt=-kkt
      isig=-isig
175      llt=isydeg(ll,n)
      if(llt) 176,200,180
176   llt=-llt
      isig=-isig
180   nklst=2
       if(iit.eq.kk.and.jjt.eq.ll.and.kkt.eq.ii.and.llt.eq.jj)
     &    yijkl=.true.
      if(iit.ge.jjt) go to 181
      nnt=iit
      nklst=nklst+1
      iit=jjt
      jjt=nnt
181   if(kkt.ge.llt) go to 182
      nnt=kkt
      kkt=llt
      nklst=nklst+2
      llt=nnt
182   if(iit-kkt) 183,184,185
183   call ex(iit,kkt,jjt,llt)
      nklst=nklst+4
      go to 185
184   if(llt.gt.jjt) go to 183
185   ijt=iit*(iit-1)/2+jjt
      if(ijt.gt.nij) go to 186
      klt=kkt*(kkt-1)/2+llt
      if(ijt.lt.nij) go to 187
      if(klt.lt.nkl) go to 187
186   ijst=its(iit,jjt)
      klst=its(kkt,llt)
       if(its(ijst,klst).ne.1) go to 500
      go to 200
187   continue
      ijklt=ijt*(ijt-1)/2+klt
       if(ktypsy(ijklt).ne.10.and.ktypsy(ijklt).ne.1) go to 222
      if(nklst.eq.5) nklst=10
      if(nklst.ne.10) nklst=1
222   continue
      ktypsy(ijkls)=ltabs(nklst,ktypsy(ijklt))
      nad(ijkls)=isig*nad(ijklt)
      if(yprtv) write(6,190) ii,jj,kk,ll,iit,jjt,kkt,llt,nij,nkl,ijkls,
     *nad(ijkls),ktypsy(ijkls)
191   format(1x,20i5)
190   format(1x,4i3,'equivalent a',4i3,'***',2i5,'ijkls,nad,ktyp',3i7)
      go to 500
200   continue
205     continue  
      if(.not.yr) go to 220
      msis=ii
c
      if(.not.ydp) then
           do 215 i=1,nijkl
215        bijkls(i)=0.
           call transs(msis,jj,kk,ll,bijkls)
         else
           do 216 i=1,nijkl
216        bijkld(i)=0.
           call transd(msis,jj,kk,ll,bijkld,xx_pqrs,ix_pqrs)
         endif            
      kksvn=kk
      llsvn=ll
      yr=.false.
      ysaut=.false.
      if(msis.eq.0) ysaut=.true.
      monowr=0
      rewind 40
c     call chrono(tod,cput,reste,indtim)
c     ntime=cput*1.d-4
      if(yprtv)write(6,64) msi,msj,msk,msl,ntime,nijkl
64    format('0',25('*'),'pqrs sym:',4i5,i10,' /100 s     npqrs=',i7)
220   continue
      if(ysaut) go to 500
      npair=npair+2    
      ipair(npair-1)=nij     
      ipair(npair)=nkl      
      ktypsy(ijkls)=2     
      if(yijkl)      ktypsy(ijkls)=1   
      nad(ijkls)=npair/2          
      if(.not.ydp)then
           call tpqkls(nbp(ll),nbp(kk),nbp(jj),nbp(ii),nbo(jj),nbo(ii),
     &     yrs,ypq,ypqrs,yij,yijkl,ind(jj),ind(ii),bijkls)
       else
           call tpqkld(nbp(ll),nbp(kk),nbp(jj),nbp(ii),nbo(jj),nbo(ii),
     &     yrs,ypq,ypqrs,yij,yijkl,ind(jj),ind(ii),bijkld)
       endif 
c      call chrono(tod,cput,reste,indtim)
c     ntime=cput*1.d-4
      npqkl=mono2*mij
65    format(' ijkl sym',4i5,i10,' /100 s **paire:',3i4,'ijkls,nad,ktyp'
     &,3i7,'  npqkl=',i7)
      if(yprtv)write(6,65) ii,jj,kk,ll,ntime,npair,nij,nkl,ijkls,nad
     &(ijkls),ktypsy(ijkls),npqkl
70    format(1x,15f8.4)
        ncpqkl=ncpqkl+1
	if(ncpqkl.gt.idinir) then
	   write(6,*)'depassement de tableau dans itijkl'
	   write(6,*)'ncpqkl > idinir  idinir=',idinir
	   stop
        endif
        ivnir(ncpqkl,1)=npqkl
        ivnir(ncpqkl,2)=nbo(ll)
        ivnir(ncpqkl,3)=nbo(kk)
        ivnir(ncpqkl,4)=nbo(jj)
        ivnir(ncpqkl,5)=nbo(ii)
        ivnir(ncpqkl,6)=ind(ll)
        ivnir(ncpqkl,7)=ind(kk)
        ivnir(ncpqkl,8)=ind(jj)
        ivnir(ncpqkl,9)=ind(ii)
        yvnir(ncpqkl,1)=ykl
        yvnir(ncpqkl,2)=yij
        yvnir(ncpqkl,3)=yijkl         
500   continue
         if(ncpqkl.eq.0) go to 600
      if(monowr.ne.0.and..not.ydp) write(40) bufs2
      if(monowr.ne.0.and.ydp) write(40) bufd2
      rewind 40
c
           inpqkl=lsize
          do 550 i=1,ncpqkl
c
          npqkl=ivnir(i,1)
      iwmij=mijkl+1
      if(.not.ydp)then
c   simple precision
         do 540 n=1,npqkl
           inpqkl=inpqkl+1
           if(inpqkl.le.lsize) go to 540
           read(40) bufs2
           inpqkl=1
540        bijkls(n)=bufs2(inpqkl)
c     write(6,*) i,(ivnir(i,jkl),jkl=1,9),(yvnir(i,jkl),jkl=1,3)
          call tijkls(ivnir(i,2),ivnir(i,3),ivnir(i,4),
     &ivnir(i,5),ivnir(i,6),ivnir(i,7),ivnir(i,8),ivnir(i,9),
     &yvnir(i,1),yvnir(i,2),yvnir(i,3),nbp(llsvn),nbp(kksvn),yrs
     &,bijkls)
         else
c  double precision
        do 545 n=1,npqkl
           inpqkl=inpqkl+1
           if(inpqkl.le.lsize) go to 545
           read(40) bufd2
           inpqkl=1
545        bijkld(n)=bufd2(inpqkl)
c     write(6,*) i,(ivnir(i,jkl),jkl=1,9),(yvnir(i,jkl),jkl=1,3)
          call tijkld(ivnir(i,2),ivnir(i,3),ivnir(i,4),
     &ivnir(i,5),ivnir(i,6),ivnir(i,7),ivnir(i,8),ivnir(i,9),
     &yvnir(i,1),yvnir(i,2),yvnir(i,3),nbp(llsvn),nbp(kksvn),yrs
     &,bijkld)
         endif
550   continue
c     call chrono(tod,cput,reste,indtim)
c     ntime=cput*1.d-4
      if(yprtv) write(6,552) ntime,mblk,mijkl
  552 format(' fin transf. pqkl en ijkl ntime,mblk,mijkl',3i10)
600   continue
      npaire=npair/2
60    format(1x,20i6)
      write(6,*)'npaire,npair',npaire,npair
      write(25) npaire,(ipair(i),i=1,npair)
      write(25) nsyt,(nad(i),i=1,nsyt)
      write(25) (ktypsy(i),i=1,nsyt)
      if(mijkl.ne.0.and..not.ydp) call wtrbls
      if(mijkl.ne.0.and.ydp) call wtrbld
c        call chrono(tod,cput,reste,indtim)
c     ntime=cput*1.d-6
c        ntim=reste
c        write(6,610) ntime,ntim
 610     format('0fin ijkl temps ecoule :',i10,'s;reste',i10,'s')
      go to 695
690   ii=nsym+1  
      if(.not.ydp)call transs(ii,nsym,nsym,nsym,bijkls)
      if(ydp)call transd(ii,nsym,nsym,nsym,bijkld,xx_pqrs,ix_pqrs)
695   continue
      deallocate(xx_pqrs,ix_pqrs)
      return
      end
