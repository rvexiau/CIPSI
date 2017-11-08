      subroutine rotat1 (nveci,nvecf,noccup,isurci,ipass,
     1 nsym1,nsym2,nsymp,yrot,ygel,s,vecsym,srec,a,csc,yf11)
      implicit real*8 (a-h,o-x,z),logical*1(y)
      include 'pshf.prm'
      common norb,nsym,its(nsymz,nsymz),itsym(doa),itsyv(2*doa),
     1 nvec,isymp(nsymz),mijkl,nrec,mblk
      common isydeg(nsymz,nsymz),ntrsy,istop,iarwr(2*doa),
     1 nocsy(nsymz),ymono,yprt,yprtv,yatom,ywvect,ynocsy
      common /om/ c(doa*doa*2),nbp(nsymz),nbo(nsymz),ind(nsymz+1),
     y ivecsy(2*doa*nsymz)
      common /hondo/ nt,ishell(doa),kloc(doa),ktyp(doa),
     1 numc(doa),newsh(doa,nsymz),iptr(3,nsymz),idtr(6,nsymz),
     2 iftr(10,nsymz),igtr(15,nsymz),yhondo,yprth
      common/nuscf/iscf(doa),iflab(3,doa)
      dimension s(*),
     * vecsym(*),a(*)
      dimension ygel(2*doa)
      dimension srec(ns3),csc(2*doa*doa),
     * isym(2*doa),as(doa),itr(nsymz),
     1 ysigi(nsymz)
      common/dim/ns1,ns2,ns3,ns4,ns5
      numnu(i)=(i*(i-1))/2
c     write(6,*)'entree rotat1'
      nbom=norb
      nnorb=norb*norb
      nnorb1=nnorb+1
c    si yf11=true et yrot=true lecture directe des om sur 11
      if(yf11.and.yrot)then
        if(ipass.eq.1)read(11)(c(i),i=1,nnorb)
        if(ipass.eq.2)read(11) (c(i+nnorb),i=1,nnorb) 
      else
c    si yf11=false et yrot=true lecture des om sur 2(defaut)
      if(.not.yf11.and.yrot) then
        if(ipass.eq.1)then
        do 1001 ijump=1,5
 1001    read(2)
        call readb(c,nnorb,2)
        else
        do 1002 ijump=1,2
 1002   read(2)
        call readb(c(nnorb+1),nnorb,2)
        end if
c    si yrot=false lecture des om d'un passage precedent sur 11
       end if
       if(yprt) then
       ifin=0
          do iii=1,norb
             ideb=ifin+1
	     ifin=ideb+norb-1
	     write(6,*)'orbitale',iii
   	  write(6,2521)(c(k),k=ideb,ifin)
          enddo
       endif
 2521  format(10(1x,f10.7))
      if(.not.yrot)then
         if(ipass.eq.1)call readb (c,nnorb,11)
         if(ipass.eq.2)call readb (c(nnorb1),nnorb,11)
      end if
      end if
      ideb=0
      if(ipass.eq.2) ideb=nnorb
      nacti=0
      do 1006 i=nveci,nvecf
c      write(6,*) 'i,gel(i)',i,ygel(i)
      if(ygel(i))go to 1006
      nacti=nacti+1
      do 1007 j=1,norb
      c(ideb+(nacti-1)*norb+j)=c(ideb+(i-1)*norb+j)
 1007 continue
 1006 continue
      nvecf=nvecf-(norb-nacti)
c      write(6,*) 'nacti',nacti
       k=0
       do 7 i=nsym1,nsym2
      do 7 j=1,norb
      k=k+1
7     vecsym(k)=ivecsy(j+(i-1)*norb)
      do 6 i=nveci,nvecf
6     isym(i)=1
      isymp(1)=1
      if(nsym1.eq.nsym2) go to 1500
      if(yatom) go to 10
      ii=norb*(norb+1)/2
      rewind 2
      read(2)
      call readb(s,ii,2)
	 write(6,*)'matrice de recouvrement'
	 ifin=0
	 do i=1,norb
	    ideb=ifin+1
	    ifin=ideb+i-1
	    write(6,'(10f8.4)')(s(k),k=ideb,ifin)
         enddo
 1015 go to 13
608   format(' matrice unite pour s')
10    k=0
      write(6,608)
           do 12 i=1,norb
        do 11 j=1,i
        k=k+1
11      s(k)=0.d0
12      s(k)=1.d0
  13  if(yprt)write(6,601)
601    format(' vecteurs de symmetrie apres  normalisation')
      do 310 i=nveci,nvecf
  310 isym(i)=itsyv(i)
      write(6,*)(isym(i),i=nveci,nvecf)
  306 k=0
      do 320 is=nsym1,nsym2
      rec=0
      kmis=(is-nsym1)*norb
      kmit=(is-1)*norb
      do 315 i=1,norb
      do 315 j=1,norb
      if(i.ge.j) ij=i*(i-1)/2+j
      if(i.lt.j) ij=j*(j-1)/2+i
315   rec=rec+vecsym(kmis+i)*s(ij)*vecsym(kmis+j)
      if(rec.lt.1.d-08) write(6,609) is
      if(rec.lt.1.d-8)rec=1.d0
609   format(' norme nulle pour le vecteur de symetrie:',i3)
      rec=1.d0/dsqrt(rec)
      do 316 i=1,norb
316   vecsym(kmis+i)=vecsym(kmis+i)*rec
       do 321 i=1,norb
      ivecsy(kmit+i)=0
      if(dabs(vecsym(kmis+i)).lt.1.d-08) go to 321
      ivecsy(kmit+i)=is
321   continue
      if(yprtv)write(6,61) (ivecsy(kmit+i),i=1,norb)
320   continue
c
      nsymp1=nsym1+1
      if(ipass.eq.2) nsymp1=nsym1
      do 150 is=nsymp1,nsym2
      kmis=(is-1)*norb
      ism=is-1
      do 120 js=1,ism
      kmjs=(js-1)*norb
      do 110 i=1,norb
      if(ivecsy(kmjs+i).eq.0.and.ivecsy(kmis+i).ne.0) go to 120
       if(ivecsy(kmis+i).eq.0.and.ivecsy(kmjs+i).ne.0) go to 120
110   continue
       isymp(is)=isymp(js)
       go to 150
120   continue
        nsymp=nsymp+1
      isymp(is)=nsymp
150   continue
      iss=0
      do 170 is=1,nsymp
      do 160 js=is,nsym
      if(isymp(js).ne.is) go to 160
      iss=iss+1
      if(js.eq.iss) go to 160
      isymp(js)=isymp(iss)
      isymp(iss)=is
      do 155 i=1,norb
      if(isym(i).ne.js) go to 154
      isym(i)=iss
      go to 155
154   if(isym(i).ne.iss) go to  155
      isym(i)=js
155   continue
      kmis=(js-1)*norb
      kmjs=(iss-1)*norb
      do 147 i=1,norb
      a(i)=vecsym(kmis+i)
      vecsym(kmis+i)=vecsym(kmjs+i)
      vecsym(kmjs+i)=a(i)
      idum=ivecsy(kmis+i)
      ivecsy(kmis+i)=ivecsy(kmjs+i)
      ivecsy(kmjs+i)=idum
 147  continue
160   continue
170   continue
180    continue
      write(6,610) nsymp
610   format(' nombre de symetries pour les orb atomiques',i4)
      if(yprt)write(6,61) (isymp(i),i=nsym1,nsym2)
      if(yprt)write(6,602)(i,i=nsym1,nsym2)
602       format('0symetrie',15(i5,3x))
      do  603 i=1,norb
      if(yprt)write(6,604)i,(vecsym(i+(j-nsym1)*norb),j=nsym1,nsym2)
603      continue
604        format(1x,i3,(4x,15f8.4))
      if(itsyv(nveci).ne.0) go to 1500
c
      do 210 i=nveci,nvecf
      ni=(i-1)*norb
      do 210 k=1,norb
      rec=0.d0
      do 211 l=1,norb
      if(k.ge.l) kl=numnu(k)+l
      if(k.lt.l) kl=numnu(l)+k
  211 rec=rec+c(ni+l)*s(kl)
  210 csc(ni+k)=rec
c
      mtr=0
      do 250 i=nveci,nvecf
         ntr=0
      ni=(i-1)*norb
      isy=0
      rec=0.
      do 230 is=nsym1,nsym2
      rec=0.
      kmis=(is-nsym1)*norb
      do 216 l=1,norb
      kmis=kmis+1
216      rec=rec+vecsym(kmis)*csc(ni+l)
      reco=rec
      rec=dabs(rec)
      if(rec.lt.1.d-5) go to 230
225   ntr=ntr+1
      srec(ntr)=rec
      itr(ntr)=is
       ysigi(ntr)=reco.lt.0.d0
230    continue
      if(yprtv) write(6,604) i,recm,(c(ni+k),k=1,norb)
      if(ntr.ne.0) go to 240
235   write(6,631) i
631    format('0 le vecteur ',i4,' n a pu etre attribue a 1 sym')
      stop
240   recm=5.d-6
      mtr=0
      do 243 l=1,ntr
      if(recm.ge.srec(l)) go to 243
               mtr=l
      recm=srec(l)
243   continue
      if(mtr.eq.0) go to 235
      ysig=ysigi(mtr)
      srec(mtr)=0.d0
      isy=itr(mtr)
      if(yprtv) write(6,*) isy,ysig,recm,ntr,(itr(l),l=1,ntr)
246   isym(i)=isy
60       format(1x,15f8.4)
      kmis=nbom*(isy-1)
      do  330 k=1,norb
      as(k)=c(ni+k)
      if(ivecsy(kmis+k).eq.0) as(k)=0.d0
330    continue
      if(i.eq.nveci) go to 337
      im=i-1
      do 335 j=nveci,im
      nj=(j-1)*norb
      rec=0.d0
      do 334 l=1,norb
334   rec=rec+csc(nj+l)*as(l)
      if(dabs(rec).gt.1.d-5) go to 240
335   continue
337   continue
      recm=0.
      do  340 k=1,norb
      do  340 l=1,norb
      if(k.ge.l) kl=numnu(k)+l
      if(k.lt.l) kl=numnu(l)+k
340   recm=recm+as(k)*as(l)*s(kl)
      if(dabs(recm).lt.2.d-5) go to 240
      recm=1.d0/dsqrt(recm)
      if(ysig) recm=-recm
      do  345 k=1,norb
345   c(ni+k)=recm*as(k)
c
      do 258 l=1,norb
      rec=0.d0
      do 257 k=1,norb
      if(k.ge.l) kl=numnu(k)+l
      if(k.lt.l) kl=numnu(l)+k
  257 rec=rec+c(ni+k)*s(kl)
  258 csc(ni+l)=rec
c
      isy=0
      recm=0.d0
      do 360 is=nsym1,nsym2
      rec=0.d0
      kmis=(is-nsym1)*norb
      do 350 l=1,norb
      kmis=kmis+1
350   rec=rec+vecsym(kmis)*csc(ni+l)
      if(dabs(rec).lt.dabs(recm)) go to 360
      recm=rec
      isy=is
360   continue
      isym(i)=isy
      if(recm.gt.0.d0) go to 249
      do 365 k=1,norb
      csc(ni+k)=-csc(ni+k)
365   c(ni+k)=-c(ni+k)
c
  249 continue
      if(yprtv) write(6,604) isy,(c(ni+k),k=1,norb)
250       continue
1500  continue
      kk=0
c
      koc=0
         kvirt=noccup
      do 100 is=nsym1,nsym2
      nok=0
      do 90 i=nveci,nvecf
      ni=(i-1)*norb
      if(isym(i).ne.is) go to 90
      if(ynocsy) go to 101
      if(i.gt.noccup) go to 80
      go to 102
  101 nok=nok+1
      if(nok.gt.nocsy(is)) go to 80
  102 continue
      koc=koc+1
      k=koc
      go to 81
80     kvirt=kvirt+1
      k=kvirt
   81 nk=(k-1)*norb
      do 85 j=1,norb
85     csc(nk+j)=c(ni+j)
      itsyv(k)=isym(i)
      kk=kk+1
      iscf(kk)=i
90     continue
100    continue
61       format(1x,40i 3)
      write (6,64) (itsyv(i),i=nveci,nvecf)
64    format(' symetrie des om',(1x,30i3))
c
      do 405 i=nveci,nvecf
      ni=(i-1)*norb
      do 405 k=1,norb
  405 c(ni+k)=csc(ni+k)
c
      return
      end
