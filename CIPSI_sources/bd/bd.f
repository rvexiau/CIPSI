      program bd
      implicit real*8(a-h,o-p,r-z),logical*4(q)
      character*26 timeday
      include 'bd.prm'
      character*20 title
      common/gest/spert,eold(netm),cold(ncpm,netm),ind(ncpm),qz(ncpm)
      integer*4 isytr
      integer*4 itab(lsize),jtab(lsize)
      common/info/ncf,ncper,nca,qnorm,qwa,qwb,mbuf,metat2,netat,niter
     *,qwvec,qhef,ntrsy,qrvec,seuil,qws,qsuis,qsuip,qener,njter,qgiv
     *,qictot,qsui04,qdcple,ityp
      common/serre/cser(nczm,nczm),va(nvsi),vb(nvsi),ites(netm)
      common/nom33/title
      common/tprint/tvec
      common/ceref/escf
      namelist/option/qnorm,netat,spert,ind,cold,mbuf,qwa,qwb,qrvec
     &  ,ncf,ncper,metat,test,niter,ezer,ites,qwvec,qhef,seuil,ntrsy
     &,qws,qsuis,qsuip,qener,njter,qgiv,qictot,qsui04,qdcple,tvec,title
     &,ityp
      call fdate(timeday)
      write(6,9999) timeday
 9999 format(/,80(1h*),/,8(10h   bd     ),/,  10x,a25,/,
     * 80(1h*))
      qdcple=.false.
      ezer=0.d0
           do1j=1,netm
           ites(j)=1
             do1 i=1,ncpm
1            cold(i,j)=0.d0
      qictot=.false.
      qnorm=.false.
	  tvec=0.02
          njter=1
      qwa=.false.
      qwvec=.false.
       qws=.false.
        qsuis=.false.
       qsuip=.true.
       qsui04=.false.
      qener=.true.
      qrvec=.false.
      qhef=.false.
       title='                    '
       ntrsy=0
      seuil=1.d0
      qwb=.false.
      qgiv=.true.
      mbuf=lsize
      spert=0.d0
      niter=50
      test=1.d-4
      metat=1
      ncper=1
        ind(1)=0
      netat=0
      ityp=0
      call openf
       read(4,end=7777,err=7777)ib,ib,ib,metat,ncper
       nbr=metat*(metat+1)/2
       read(4,end=7777,err=7777)(bb,i=1,nbr),escf
7777  continue
      write(6,*)
      write(6,*) ' energie de reference '
      write(6,*)escf
c     determination directe de ncf par relecture de la file 1
c     jp daudey oct 88
      rewind 1
      ncf=0
125   read(1) itab,jtab
      last=lsize
      if(itab(lsize).eq.-1)last=jtab(lsize)
      do 127 i=1,last
      if(itab(i).gt.ncf) ncf=itab(i)
      if(jtab(i).gt.ncf) ncf=jtab(i)
  127 continue
      if(itab(lsize).ne.-1) go to 125
      write(6,979)
979   format(/,100('*'),/,10('    diag'),/,120('*'),/)
      write(6,*)
      write(6,*)' methode variante de davidson.toulouse 84..pelissier'
      write(6,*)
      write(6,*)'dimension de la matrice d ic sur la file 1',ncf
      read(5,option) 
         if(ind(1).eq.0)then
         do5i=1,ncper
5        ind(i)=i
c         do10i=1,ncf
c10        scor(i)=1.d0
          write(6,980) ncf,ncper,metat 
980   format(' dimension de la matrice de moyen: ',i4,/,
     *' dimension de la matrice de depart: ',i4,/,
     *' nombre de vecteurs propres desires: ',i4)
        metat2=metat
        if(ncper.gt.nczm)then
        write(6,*)' matrice de depart trop grande'
        write(6,*)' diminuer ncper'
        stop
        endif
        if(metat.gt.netm)then
        write(6,*)' trop de vecteurs propres'
        write(6,*)' diminuer metat'
        stop
        endif
          write(6,*)' ites',(ites(i),i=1,metat)
          write(6,982) test*10.d0
982     format(1x,'precision',f10.6,/)
          endif
          do20 i=1,metat
          if(ites(i).eq.1)then
        ites(i)=0.1d0/test
         netat=netat+1
          endif
20        continue
      if(ncper.gt.nczm)then
        write(6,*)' matrice de depart trop grande'
      write(6,*)' ncper trop grand'
      stop
      endif
      iwant=ncf*netat
       iwant=max(iwant+iwant,ncper*ncper+iwant)
      write(6,981)
981   format(/,1x,100('*'),/)
        call storic(iwant)
        
      call fdate(timeday)
      write(6,9998) timeday
 9998 format(/,80(1h*),/,8(10h fin bd   ),/, 10x,a25,/,
     * 80(1h*))
      stop
      end
