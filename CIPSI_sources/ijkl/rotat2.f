      subroutine rotat2 (nveci,nvecf,noccup,isurci,ipass, 
     1 nsym1,nsym2,nsymp,yrot,ygel,inpv,trec,tdeg,s,
     2 vecsym,srec,a,ch,csc,yf11)
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
     *vecsym(*),a(*),ch(*)
      dimension ygel(2*doa)
      dimension isym(2*doa),csc(2*doa*doa)
      dimension ideg(doa),iniv(doa)
      dimension srec(ns3),itr(nsymz)
      common/dim/ns1,ns2,ns3,ns4,ns5
      numnu(i)=(i*(i-1))/2
      nbom=norb
      nnorb=norb*norb
      nnorb1=nnorb+1
c    si yf11=true et yoot=true lecture directe des om sur 11
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
c      nvecf=nvecf-(norb-nacti)
c      write(6,*) 'nacti,nvecf',nacti,nvecf
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
c  determination de la degenerescence
c
      do 700 i=1,nnorb
700   csc(i)=c(i)
      yfi50=.false.
c     write(6,*) 'avant tmono ' ,nvecf
      call tmono(ipass,nveci,nvecf,ch,csc,yfi50)
      i=1
      niv=0
710   ndeg=1
      if(i.gt.nvecf) go to 730
      niv=niv+1
      j=i
      jj=i*(i+1)/2
      cii=ch(jj)
720   j=j+1
      jj=jj+j
      if(j.gt.nvecf) then
      ideg(niv)=ndeg
      iniv(niv)=i
      go to 730
      endif
      deps=cii-ch(jj)
      deps=dabs(deps)
      if(deps.lt.tdeg) then
      ndeg=ndeg+1
      go to 720
      else
      ideg(niv)=ndeg
      iniv(niv)=i
      i=i+ndeg
      go to 710
      endif
730   continue
      write(6,*) ' nombre de niveaux ',niv
      write(6,*) ' premiere orbitale '
      write(6,755) (iniv(i),i=1,niv)
      write(6,*) ' degenerescence '
      write(6,755) (ideg(i),i=1,niv)
755   format(20(1x,i3))
      do 780 i=1,nnorb
780   c(i)=csc(i)
c     goto 781
c
c  rotation des om
c
      do 210 i=nveci,nvecf
      ni=(i-1)*norb
      if(yprt) then
      write(6,*) i
      write(6,221) (c(ni+ii),ii=1,norb)
      endif
      do 210 k=1,norb
      rec=0.d0
      do 211 l=1,norb
      if(k.ge.l) kl=numnu(k)+l
      if(k.lt.l) kl=numnu(l)+k
      rec=rec+c(ni+l)*s(kl)
c     if(i.eq.1)
c    * write(6,*)'rec',k,l,kl,numnu(k),numnu(l),ni,c(ni+l),s(kl)
  211 continue
  210 csc(ni+k)=rec
c 781 continue
c
      do 225 ni=1,niv
      i=iniv(ni)
      mi=(i-2)*norb
      ndeg=ideg(ni)
      mis=-norb
      mir=-ndeg
      ntr=0
      do 235 is=nsym1,nsym2
      mis=mis+norb
      mir=mir+ndeg
      mj=mi
      yrec=.false.
      do 230 j=1,ndeg
      mj=mj+norb
      rec=0.
      do 216 l=1,norb
c	 write(6,*)'recou',mis+l,mj+l,vecsym(mis+l),csc(mj+l)
216      rec=rec+vecsym(mis+l)*csc(mj+l)
      srec(mir+j)=rec
      rec0=dabs(rec)
      if(rec0.ge.trec) yrec=.true.
230    continue
      if(yrec) then
      ntr=ntr+1
      itr(ntr)=is
      endif
235    continue
      if(ntr.gt.ndeg.or.(ntr.lt.ndeg.and.ntr.ne.1))then
      write(6,*)' probleme de degenerescence orbitale ',i,
     *'ndeg=',ndeg,'recouvrement avec ',ntr,' vecteurs de symetrie',
     *(itr(ii),ii=1,ntr)
      write(6,*) ' orbitales'
      write(6,*)'ntr,ndeg',ntr,ndeg
      i1=mi+1
      do 218 jj=i,i+ndeg-1
      i1=i1+norb
      i2=i1+norb-1
      write(6,*) jj
218   write(6,221) (c(ii),ii=i1,i2)
      write(6,*) ' recouvrement'
      mi=-ndeg+1
      do 219 j=nsym1,nsym2
      mi=mi+ndeg
      mf=mi+ndeg-1
219   write(6,221) (srec(jj),jj=mi,mf)
 221  format(10(1x,f25.17))
      write(6,'(a17,f15.7)') ' !! ERROR in ijkl',trec
      if(ntr.lt.ndeg) then
        write(6,*) 'decrease trec ?, error tag ='
        write(6,'(i2)') 2
        stop
      else
        write(6,*) 'fix ? increase trec ?, error tag ='
        write(6,'(i2)') 1
        stop
      endif
      else
      if(ntr.lt.ndeg.and.ntr.eq.1) then
      write(6,3246) (jj,jj=i,i+ndeg-1)
 3246 format(' attention les om suivantes sont degenerees
     * accidentellement et ne sont pas tournees',/,
     * ' verifier que c''est bien ce que vous souhaitez',/,20i4)
      do 3247 ntr=1,ndeg
      isym(i+ntr-1)=itr(1)
 3247 continue
      else
c
c  projection
c
      lk=0
      do 800 ntr=1,ndeg
      is=itr(ntr)
      mir=(is-1)*ndeg
      do 810 l=1,norb
      lk=lk+1
      aa=0.d0
      jl=mi+l
      do 855 j=1,ndeg
      jl=jl+norb
 855  aa=aa+srec(mir+j)*c(jl)
810   a(lk)=aa
800   isym(i+ntr-1)=is
      do 825 i=1,ndeg
      ici=mi+norb
      iaf=ndeg*norb
      do 825 k=1,iaf
825   c(ici+k)=a(k)
      endif
      endif
225   continue
c
c    verification de l orthonormalite
c
      do 910 i=nveci,nvecf
      ni=(i-1)*norb
      do 910 k=1,norb
      rec=0.d0
      do 911 l=1,norb
      if(k.ge.l) kl=numnu(k)+l
      if(k.lt.l) kl=numnu(l)+k
  911 rec=rec+c(ni+l)*s(kl)
  910 csc(ni+k)=rec
c
c   normalisation
c
      mi=-norb
      do  345 i=nveci,nvecf
      mi=mi+norb
      rec=0.d0
      do 350 k=1,norb
350   rec=rec+c(mi+k)*csc(mi+k)
      rec0=dabs(rec)
      rec0=dsqrt(rec0)
      if(rec0.lt.1.d-6) then
      write(6,*) 'norme trop petite pour le vecteur i',i,
     *' norme=',rec0
      write(6,'(a17,f15.7)') ' !! ERROR in ijkl',trec
      write(6,*) 'fix ?, error tag ='
      write(6,'(i2)') -1
      stop
      endif
      rec0=1.d0/rec0
      do  345 k=1,norb
      csc(mi+k)=csc(mi+k)*rec0
345   c(mi+k)=rec0*c(mi+k)
c
c  verification orthogonalisation
c
      mi=-norb
      do 942 i=nveci,nvecf
      mi=mi+norb
      mj=-norb
      if(nveci.lt.i) then
      do 940 j=nveci,i-1
      mj=mj+norb
      rec=0.d0
      do 950 k=1,norb
950   rec=rec+c(mj+k)*csc(mi+k)
      rec0=dabs(rec)
      if(rec0.gt.trec) then
      write(6,*) 'les vecteurs ',i,' et',j,'ne sont pas orthogonaux',
     *' rec=',rec0
      write(6,'(a17,f15.7)') ' !! ERROR in ijkl',trec
      write(6,*) 'increase trec ?, error tag :'
      write(6,'(i2)') 1
      stop
      endif
940   continue
      end if
  942 continue
1500  continue
c
      kk=0
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
c     write(6,*) (c(i),i=1,nnorb)
      return
      end
