!
!     DB (06/2014) : 
!       Appels à 'giveis' remplacés par des appels à 'diagonaliser'
!

      subroutine  model(dtab)
      implicitreal*8(a-h,o-p,r-z),logical*4(q)
      integer smax      
      include 'bd.prm'
      logical*1 qion
      integer*4 nepp,trpp
      common/gest/spert,eold(netm),cold(ncpm,netm),ind(ncpm),qz(ncpm)
      common/info/ncf,ncper,nca,qnorm,qwa,qwb,mbuf,metat2,netat,niter
     *,qwvec,qhef,ntrsy,qrvec,seuil,qws,qsuis,qsuip,qener,njter,qgiv
     *,qictot,qsui04,qdcple
      common/serre/cser(nczm,nczm),va(nvsi),vb(nvsi),ites(netm)
      dimension dtab(1),cper(ncpm,netm),c(ncpm*9),qsel(nczm)
      dimension numero(ncpm)
       common/wry/ work(ncpm*(ncpm+1)/2),cw(ncpm*ncpm),ew(ncpm)
      equivalence(cper(1,1),work(1)),(c(1),va(1))
      
c      
c      Remplissage de la moitié de matrice
c
      ii=0
      do 420 i = 1, nca
          do 410 j = 1, i
              ii=ii+1
410           work(ii)=dtab(ii)
420   continue

      nd=nca
      
c               Si premier appel depuis bigdia :      
                if(nca.eq.ncper.and.netat.eq.0)then
      if(qsui04) then
      rewind 4
      write(6,*)' qsui04=t :  vecteurs de depart lus sur la file 04 '
      read(4,err=998) nomr,noca,nocb,nvr,ndetr,isz,qion,ii,ik,
     *(nepp,j=1,2*ndetr),(trpp,j=1,2*ii),
     *(cw(j),j=1,ik),(ew(j),j=1,nvr)
      goto 997
 998  write(6,*) ' !!erreur dans la lecture de la file 04 dans model!!'
 997  continue
      write(6,*) ' '
      write(6,*) ' metat04=',nvr,'      metat04*ncf04=',ik
      else
         call diagonaliser( nd, work, ew,cw ,nd )
         !call giveis(nd,nd,nd,work,c,ew,cw,ierr)
      endif
          mnc=0
          do407j=1,nca
          do406i=1,j
             mnc=mnc+1
406       work(mnc)=dtab(mnc)
407       work(mnc)=dtab(mnc)*.5d0
       do454i=1,ncper
454      qsel(i)=.false.
      do460j=1,metat2
       if(ites(j).le.0)goto460
        netat=netat+1
        ites(netat)=ites(j)
           eold(netat)=ew(j)
       t=seuil
      do455i=1,ncper
      cser(i,netat)=cw(i+ncper*(j-1))
455    t=max(t,abs(cser(i,netat)))
       t=t*seuil
        do456i=1,ncper
456     if(abs(cser(i,netat)).gt.t)qsel(i)=.true.
460      continue
       nd=netat
      do480j=1,metat2
       if(ites(j).ne.0)goto480
        netat=netat+1
        ites(netat)=ites(j)
           eold(netat)=ew(j)
      do475i=1,ncper
475   cser(i,netat)=cw(i+ncper*(j-1))
480      continue
        metat2=netat
        netat=nd
         mnc=0
         do205nc=1,metat2
         do204mc=1,nc
          mnc=mnc+1
204      dtab(mnc)=0.d0
205      dtab(mnc)=eold(nc)
        do550i=1,ncper
        if(.not.qsel(i))goto550
        rec=1.d0
        do510j=1,metat2
510     rec=rec-cser(i,j)*cser(i,j)
        if(rec.lt.seuil)goto550
         if(metat2.eq.netm)goto550
         metat2=metat2+1
        rec=1.d0/sqrt(rec)
        do530j=1,ncper
        t=0.d0
        do520k=1,metat2-1
520     t=t-cser(j,k)*cser(i,k)*rec
        cser(j,metat2)=t
530     continue
          cser(i,metat2)=cser(i,metat2)+rec
        write(6,*)' supplement:',metat2
        write(6,100)(cser(j,metat2),j=1,ncper)
                  do534j=1,ncper
534               vb(j)=0.d0
                 jk=0
                 do535j=1,ncper
                 do535k=1,j
                 jk=jk+1
                 vb(k)=vb(k)+work(jk)*cser(j,metat2)
535              vb(j)=vb(j)+work(jk)*cser(k,metat2)
         do545j=1,metat2
          mnc=mnc+1
         t=0.d0
         do540l=1,ncper
540       t=t+cser(l,j)*vb(l)
545       dtab(mnc)=t
550      continue
c                do405i=1,nca
c                if(cold(i,1).ne.0.d0)goto407
c405             continue
      write(6,*) ' diagonalisation de la matrice de depart'
      write(6,*)
      write(6,*) ' valeurs propres'
      write(6,100)(ew(i),i=1,metat2)
                   do197k=1,metat2
                    do196i=1,ncper
196         cold(i,k)=0.d0
197          cold(k,k)=1.d0
               mnc=0
           return
c     Si ce n'est pas le premier appel depuis bigdia :           
      else
         call diagonaliser( nd, work, ew, cw,nd )
         !call giveis(nd,nd,nd,work,c,ew,cw,ierr)
            endif
c           Fin du "if" premier appel ou non

c407         continue
      do430k=1,netat
      ji=0
      do425j=1,nca
       rec=0.d0
       if(qener)then
        if(k.eq.j)rec=1.d0
       elseif(qsuis)then
             rec=dabs(cw(ji+k))
             ji=ji+nca
           else
       do424i=1,nca
      ji=ji+1
      rec=rec+cold(i,k)*cw(ji)
424       continue
        endif
        cper(j,k)=abs(rec)
425       continue
430       continue
      do440k=1,netat
434      recm=0.d0
         recn=0.d0
      do435j=1,nca
      if(recm.lt.cper(j,k))then
      recn=recm
      recm=cper(j,k)
       jjs=js
      js=j
      else if(recn.lt.cper(j,k))then
       recn=cper(j,k)
       jjs=j
      endif
435         continue
       if(recm.lt.0.5d0)then
       write(6,*)' le vecteur ',k,' est douteux rec ',recm
               goto438
      endif
        do436i=k+1,netat
        if(cper(js,i).gt.recm)then
           recp=0.d0
           do 433j=1,nca
        if(cper(j,i).gt.recp.and.j.ne.js)then
          recp=cper(j,i)
          endif
433     continue
              if(recp.lt.recn)then
         js=jjs
          recn=0.d0
        endif
         endif
436       continue
438        continue
          do437i=1,nca
437       cold(i,k)=cw(i+nca*(js-1))
c           metat2=max(js,metat2)
           eold(k)=ew(js)
          do439i=1,netat
439        cper(js,i)=0.d0
440       continue
      ii=0
      i1=1
      return
      entry final
      nn=0
      maxi=0
29    maxi1=maxi+1
      maxi=maxi1+9
      if(metat2.lt.maxi) maxi=metat2
      write(6,1006)
      do 12 k=1,nca
      write(6,100) (cw(nca*(j-1)+k),j=maxi1,maxi)
12      continue
100   format(1x,10f12.6)
      if(maxi.lt.metat2) go to 29
1006  format(1x,'vecteurs propres')
      return
      entry   prt(dtab,n1,n2)
       do800i=n1,n2
800      write(6,801)(dtab(i*(i-1)/2+j),j=1,i)
801     format(1x,12f10.5)
       return
      end
