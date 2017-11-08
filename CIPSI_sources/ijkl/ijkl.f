      program ijkl
      implicit real*8 (a-h,o-x,z),logical*1(y)
      include 'pshf.prm'
      character *26 timeday  
      real*4 bijkl
      character*8 group,grp1,cinfv_min,dinfh_min,cinfv_maj,dinfh_maj,
     * bname,typpro
      logical chstan
        dimension ardp(kget/2)
c     le tableau ardp est exactement identique a bijkl
c     mais declare en double precision
c
      dimension ygel(2*doa)
      dimension title(10),ioda(19)
      common norb,nsym,its(nsymz,nsymz),itsym(doa),itsyv(2*doa),
     1 nvec,isymp(nsymz),mijkl,nrec,mblk
      common isydeg(nsymz,nsymz),ntrsy,istop,iarwr(2*doa),
     1 nocsy(nsymz),ymono,yprt,yprtv,yatom,ywvect,ynocsy
      common /om/ c(2*doa*doa),nbp(nsymz),nbo(nsymz),ind(nsymz+1),
     y ivecsy(2*doa*nsymz)
      common /hondo/ nt,ishell(doa),kloc(doa),ktyp(doa),
     1 numc(doa),newsh(doa,nsymz),iptr(3,nsymz),idtr(6,nsymz),
     2 iftr(10,nsymz),igtr(15,nsymz),yhondo,yprth
      common/rot/csc(2*doa*doa)
      common/dim/ns1,ns2,ns3,ns4,ns5
      common/bij/bijkl(kget)
        equivalence (ardp(1),bijkl(1))
      namelist/ijklin/nsym,norb,itsym,itsyv,isymp,nvec,yhondo
     *,noc,ntrans,yprtv,istop,group,yauto,ljkf,nsyma,nsymb
     * ,nrec,mblk,mijkl,naxis,yrot,yuhf
     *,ymono,yprt,nveci,nvecf,ivecsy,yatom
     &,iwant,ycv,ywvect,yprth,ynocsy,nocsy,yorbf
     &,inpv,trec,tdeg,yrecsp,ygel,yf11,ydp,chstan,yreduc
      data cinfv_min,dinfh_min / 'cinfv','dinfh'/
      data cinfv_maj,dinfh_maj / 'CINFV','DINFH'/

      nsym=1
      nfi11=11
      nfi10=10
      call openf  
      call fdate(timeday)
      write(6,9999)timeday
      do 2 i=1,doa*nsymz
2     ivecsy(i)=0
      do 8 i=1,nsymz
    8 nocsy(i)=0
      ntrans=-1
      ntran=0
c    attention  si on impose ntrans=0 dans la namelist ,cette valeur
c    sera imposee meme si yauto=.true.
c    modif pour la version actuelle de mcscf qui ne tient pas compte
c    des equivalences x/y
c    dans ce cas uitiliser les groupes cinfv ou dinfh de maniere standar
c    pour faire la reconnaissance mais mettre ntrans=0 dans la namelist
c    attention encore au cas des donnes non automatiques pour les groupe
c    qui n'ont pas d'equivalence .specifier ntrans=0
      nsymp=1
      istop=0
      yprt=.false.
      chstan=.true.
      trec=1.d-5
      tdeg=1.d-6
      inpv=0
      yprtv=.false.
      ywvect=.false.
      ynocsy=.false.
      yrot=.true.
      ydp=.false.
      norb=0
      noc=0
      ymono=.false.
      yrecsp=.true.
      ycv=.false.
      yorbf=.true.
      group=' '
      iwant=0
      nveci=1
      nvecf=0
      nvec=0
      ljkf=0
      yhondo=.true.
      yreduc=.true.
      nsyma=0
      nsymb=0
      yprth=.false.
      yatom=.false.
      nrec=0
      mblk=0
      mijkl=0
      do 27 i=1,2*doa
      ygel(i)=.false.
   27 itsyv(i)=0
      do 5 i=1,nsymz
5     isymp(i)=i
      yf11=.false.
      yauto=.true.
      naxis=2
      read(5,ijklin,end=998)
      rewind 2
      read(2) ioda,nt,nshell,num,nveca,nvecb,bname,ngauss,ne,
     1 ich,mul,na,nb,nat,group,naxis,title
      if(group.eq.' ') group=grp1
      call rondof(norb1,noc1,enuc,yauto,group,ivecsy,
     * nsym,ntran,naxis,nshell,nat,num,ne,yf)
      if(ntrans.eq.(-1))ntrans=ntran
      if(nvecf.eq.0) nvecf=nveca
      write(6,*) 'num,nvecf,ne',num,nvecf,ne
      if(.not.yuhf) bname='rhf'
      if(bname.ne.'uhf')nvecb=0
      if(nsyma.eq.0) nsyma=nsym
    3 continue
c      if(.not.yrot) read(11)
c      if(.not.yrot) read(11)
      if(norb.eq.0) norb=norb1
      nnorb=norb*norb
      nnorb1=nnorb+1
      typpro=' pshf'
      if(norb.le.doa) go to 6
      write(6,9992) norb,doa
      stop
    6 continue
      if(noc.eq.0) noc=noc1
      if(.not.ynocsy) go to 45
      noc=0
      do 44 i=1,nsym
   44 noc=noc+nocsy(i)
   45 continue
      if(nvecf.eq.0) nvecf=norb
      yprts=yprtv
      if(.not.yprth) go to 7
      write(6,9994) (i,i=1,nsym)
      do 4 i=1,norb
      n1=i
      n2=i+(nsym-1)*norb
    4 write(6,9993) i,(ivecsy(k),k=n1,n2,norb)
    7 continue
      ns=nsym
      if(.not.yauto.and.istop.ne.1)
     *read(5,*,end=99,err=99)((its(i,j),j=1,ns),i=1,ns)
  13  format(1x,20i6)
      ntrsy=0
      if(ntrans.eq.0) go to 16
 168  format('transformation de symetrie',/,1x,20i6)
      if(yauto) j=1 
      do 15 i=1,ntrans
      if(.not.yauto.and.istop.ne.1)
     *read(5,*) j,((isydeg(k,l+ntrsy),k=1,nsym),l=1,j)
  15   ntrsy=ntrsy+j
16    continue
      ns1=norb*(norb+1)/2
      ns2=ns1+nsym*norb
      ns3=ns2+nsym*norb
      ns4=ns3+nsym*norb
      ns5=ns4+(norb*(norb+1)/2)
      iwant=ns5+nnorb
      if(yuhf) iwant=iwant+nnorb
      iget=0
c     dans la partie de reconnaissance et de reclassement
c     le tableau ardp contient
c     de iget a ns1+iget la matrice s
c     de iget+ns1 a iget+ns2 la matrice des vecteurs de symetrie
c    de iget+ns2 a iget + ns3 puis de iget+ns3 a iget+ns4 puis
c    de iget+ns4 a iget+ns5 des
c    matrices de travail
c     de ns4+iget a ns4+nnorb+iget la matrice des coefficients
c
      ipass=1
      nveci=1
c      nvecf=nveca
      nsym1=1
      nsym2=nsyma
      ns3=ns1+ns2
      isurci=1
1613  isurcj=isurci
      noccup=nveci-1      
      if(yrecsp)then
      call rotat2(nveci,nvecf,noc,isurci,ipass,
     1 nsym1,nsym2,nsymp,yrot,ygel,inpv,trec,tdeg,
     2 ardp(iget+1),ardp(iget+ns1+1),ardp(iget+ns2+1),
     3 ardp(iget+ns3+1),ardp(iget+ns4+1),ardp(iget+ns5+1),yf11)
c
      else
      call rotat1 (nveci,nvecf,noc,isurci,ipass,nsym1,nsym2,nsymp,
     1 yrot,ygel,ardp(iget+1),ardp(iget+ns1+1),ardp(iget+ns2+1),
     2    ardp(iget+ns3+1),ardp(iget+ns5+1),yf11)
      end if       
      nvec=nvec+nvecf-nveci+1
      write(6,9997) norb,nvec,nsym,ntrans,nveci,nvecf,yauto,ljkf,enuc
c    si yreduc=t, passage de 14 a 10 rep. pour dinfh
c                 passage de 7 a 5   rep. pour cinfv
      if(yreduc.and.yf.and.yorbf) then
      if(group.eq.dinfh_min.or.group.eq.dinfh_maj) then
      nsym=10
      do 8010 i=11,14
      do 8010 j=1,norb
      if(ivecsy((i-1)*norb+j).ne.0) then
      ivecsy((i-9)*norb+j)=ivecsy((i-1)*norb+j)
      end if
 8010 continue
      do 8020 i=nveci,nvecf
      if(itsyv(i).gt.10) itsyv(i)=itsyv(i)-8
 8020 continue
      nsymp=nsymp-2
      end if
      if(group.eq.cinfv_min.or.group.eq.cinfv_maj) then
      nsym=5
      do 8030 i=6,7
      do 8030 j=1,norb
      if(ivecsy((i-1)*norb+j).ne.0) ivecsy((i-5)*norb+j)=
     * ivecsy((i-1)*norb+j)
 8030 continue
      do 8040 i=nveci,nvecf
      if(itsyv(i).gt.5) itsyv(i)=itsyv(i)-4
 8040 continue
      nsymp=nsymp-2
      end if      
      write(6,8050) nsym,nsymp,(isymp(i),i=nsym1,nsym2)
 8050 format(' ****************reduction du nombre de symetries',/,
     * ' nombre de representations irreductibles  ',i3,/,
     * ' nombre de symetries pour les primitives  ',i3,/,
     * ' type de symetrie des primitives pour les r.i. ',20i3)
      end if
      if(ipass.eq.1)
     1 call wrtdb(c,nnorb,nfi10)
      if(ipass.eq.2)
     1 call wrtdb(c(nnorb1),nnorb,nfi10)
      if(ipass.eq.2.or.bname.ne.'uhf') go to 1614
      ipass=2
      nveci=norb+1
      nvecf=norb+nvecb
      nsym1=nsyma+1
      nsym2=nsyma+nsymb
      ns3=ns3+nnorb
      go to 1613
 1614 continue
c     si bname.ne.uhf on ecrit deux fois la matrice des vecteurs 
c     alfa pour garantir la compatibilite avec le reste de la chaine
      if(bname.ne.'uhf') call wrtdb (c(1),nnorb,nfi10)
c   
      call wrtdb(ardp(iget+1),norb*(norb+1)/2,nfi10)
      if(istop.ne.1) go to 1615
      stop
 1615 continue
      write(6,161)
161   format('0      table de multiplication du groupe considere')
      do 12 i=1,nsym
12    write(6,13)(its(i,j),j=1,nsym)
c
      if(yhondo.and..not.ymono)
     1           write(25) nsym,nvec,noc,ntrans,(itsyv(i),i=1,nvec),
     2 norb,ydp,chstan
      if(yhondo.and..not.ymono)
     1           write(25) ((its(i,j),i=1,nsym),j=1,nsym)
      if(ntrans.eq.0.or..not.yhondo) go to 35
      j=1
      write(6,121)
121   format('0      transformations de symetrie')
      do 34 i=1,ntrsy
      write(6,13) (isydeg(k,i),k=1,nsym)
      if(.not.ymono) write(25) j,(isydeg(k,i),k=1,nsym)
   34 continue
 35   continue 
      call vector(ardp(ns5+iget+1))         
      n=0
      m=0
      iwant=0
      do 7546 isym=1,nsym
      if(nbp(isym).le.m) go to 7546
      do 7545 jsym=1,isym
      ijs=its(isym,jsym)
      do 7544 ksym=1,isym
      lmax=ksym
      if(ksym.eq.isym) lmax=jsym
      do 7543 lsym=1,lmax
      kls=its(ksym,lsym)
      if(its(ijs,kls).ne.1) go to 7543
      nn=nbp(isym)*nbp(jsym)
      if(isym.eq.jsym) nn=(nn+nbp(isym))/2
      mm=nbp(ksym)*nbp(lsym)
      if(ksym.eq.lsym) mm=(mm+nbp(ksym))/2
      mmnn=mm*nn
      if(isym.eq.ksym.and.jsym.eq.lsym) mmnn=(mmnn+nn)/2
      no=nbo(isym)*nbo(jsym)
      if(isym.eq.jsym) no=(no+nbo(isym))/2
      mo=nbo(ksym)*nbo(lsym)
      if(ksym.eq.lsym) mo=(mo+nbo(ksym))/2
      mno=mo*no
      if(isym.eq.jsym.and.ksym.eq.lsym) mno=(mno+no)/2
      iwant=max0(iwant,mmnn,mno,nn*mo,mm*no)
 7543 continue
 7544 continue
 7545 continue
 7546 continue
10    continue
      if(ymono) iwant=1000
      write(6,910) iwant
910   format('0 taille maximum pour le bloc pqrs ou pqkl (variable)',i8)
      write(6,9996) delt
      if(istop.eq.1) stop
      if(istop.eq.2) yprtv=yprts
      yprth=yprtv
c
      if(istop.eq.2) stop
      if(ydp) then
           if(kget.lt.iwant*2) then
           write(6,*) 'la dimension kget ',kget,' de bijkl doit etre 
     *     superieure a deux fois iwant( cas double precision) ',2*iwant
           stop
           end if
       else
           if(kget.lt.iwant) then
           write(6,*) 'la dimension kget ',kget,' de bijkl doit etre 
     *     superieure a iwant( cas simple precision) ',iwant
           stop
           end if
       endif
      if(istop.eq.3) yprtv=yprts
      if(.not.ymono)call itijkl(bijkl,bijkl,ydp)
      if(yhondo.and..not.ymono) write(25) enuc,ljkf
      if(istop.eq.3) stop
      if(ymono) yprtv=yprts
      iwant=ns1+nnorb
      if(yuhf) iwant=iwant+nnorb
      ipass=1
      nveci=1
      nvecf=nveca
 1732 continue
      yfi50=.true. 
      call tmono(ipass,nveci,nvecf,ardp(iget+1),ardp(iget+ns1+1),
     1 yfi50)     
      if(ipass.eq.2.or.bname.ne.'uhf') go to 1734
      ipass=2
      nveci=norb+1
      nvecf=norb+nvecb
      go to 1732
 1734 continue

      call fdate(timeday)
      write(6,9998)timeday
 9999 format(/,80(1h*),/,8(10h   ijkl   ),/, 10x,a25,/,
     * 80(1h*))
 9998 format(/,80(1h*),/,8(10h fin ijkl ),/, 10x,a25,/,
     * 80(1h*))
      stop
  99  continue
 998  write(6,988)
 999  write(6,989)
 988  format(' erreur fin de fichier sur  namelist ijklin')
 989  format(' erreur                    sur  namelist ijklin')
 9992 format(/35h le nombre d'orbitales moleculaires,i4,72h depasse la l
     yimite actuelle du programme',i4,' ... le programme s'arrete)
 9997 format(//27h nombre d'orbitales de base,i5,/,
     y 32h nombre d'orbitales moleculaires,i5,/,
     y 40h nombre de representations irreductibles,i5,/,
     y 22h nombre de generateurs,i5,/,
     y 37h indice de la premiere om transformee,i5,/,
     y 37h indice de la derniere om transformee,i5,/,
     y 35h generation automatique des donnees,l5,/,
     y 24h impression des j,k et f,i5,/,
     y 44h energie nucleaire (+ energie des om gelees),f20.8)
 9996 format(/' donnees completes ....   temps=',f15.5)
 9995 format(/' fin du calcul  .... temps=',f15.5)
 9994 format(/' vecteurs prototypes de chaque symmtrie (generes ou lus
     &dans le tableau ivecsy)',/,5x,20i5)
 9993 format(21i5)
      stop
      end
