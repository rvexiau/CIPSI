      subroutine  bigdia(cv,hcv)
      implicitreal*8(a-h,o-p,r-z),logical*4(q)
      include 'bd.prm'
      character*40 typener
      character*20 title
      common/nom33/title
      common/gest/spert,eold(netm),cold(ncpm,netm),ind(ncpm),qz(ncpm)
      common/info/ncf,ncper,nca,qnorm,qwa,qwb,mbuf,metat2,netat,niter
     *,qwvec,qhef,ntrsy,qrvec,seuil,qws,qsuis,qsuip,qener,njter,qgiv
     *,qictot,qsui04,qdcple,ityp
      common/serre/cser(nczm,nczm),va(nvsi),vb(nvsi),ites(netm)
       common/wry/ work(ncpm*(ncpm+1)/2),cw(ncpm*ncpm),ew(ncpm)
      common/ceref/escf
      real*8,dimension(:),allocatable :: hmat
      dimension h(ncpm*(ncpm+1)/2),cv(1)
      dimension hcv(1),iperf(netm)
      logical*1 qef
      dimension qef(kgetq)
      write(6,*) 'bigdia appel getmat'
      call getmat(h,hcv,cv)
      write(*,*) '***************************************************'
      do iii=1,10
        write(*,'(20f20.12)') (h(ij),ij=iii*(iii-1)/2+1,iii*(iii+1)/2)
      enddo
!      open (unit=852,file='mattest.dat') 
!      do iii=1,ncper
!      write(852,'(1500f20.12)') (h(ij),ij=iii*(iii-1)/2+1,iii*(iii+1)/2)
!      enddo      
!      close(852)
      write(*,*) '***************************************************'  
      
c      RV 2015 : Lapack diagonalization if number of determinant is low
      if(ncper.ge.(ncf-metat2)) then
        write(6,*) '** ncf too low, Lapack diagonalization used '
        netat=metat2
        allocate(hmat(ncf*(ncf+1)/2))
        ij=1
        do i=1,ncf
          do j=i,ncf
            hmat(ij)=h(j*(j-1)/2+i)
            ij=ij+1
          enddo
        enddo  
        call diagonaliser(ncf,h,eold(1:netat),cv(1:ncf*netat),netat)
        deallocate(hmat)
        cw=0d0
        cw(1:ncf*netat)=cv(1:ncf*netat)
        write(6,6004)(eold(i),i=1,netat)
        typener='energies var.'//title
        call stkener(typener,eold,netat)
        if(qhef.or.qwvec.or.qws)then
        nca=min(ncf,ncper)
        write(6,987)
987     format(/,18('*'),/,' resultats finaux ',/,18('*'),/)
        write(6,*) ' ***  norme des etats ***'
        write(6,*) ' etat       norme dans S      norme total'       
        do mm=1,netat
          dnorme=0d0
          dnormes=0d0
          do i=1,nca
            dnorme= dnorme + cv((mm-1)*ncf+i)*cv((mm-1)*ncf+i)
          enddo       
          dnormes=dnorme
          write(6,'(1x,i3,4x,f15.8,2x,f15.8)') mm,dnormes,dnorme
        enddo  
        write(6,*)
          write(6,*)' composantes des vecteurs propres dans s'
          if(qwb)call final
      if(qws)write(55)ncper,netat,(cv(k),k=1,ncf*netat),
     *(eold(k),k=1,netat)
           rewind 59
           write(59) ncf,netat,(cv(i),i=1,ncf*netat),
     *            (eold(i)+escf,i=1,netat)
           if(qhef) then
           if(ityp.ne.0) then 
             do i =1,3
                  write(6,*) 'bigdia  appel qhef',ityp,i,ncf,netat
              call extr(ncf,netat,cv,hcv,eold,qictot,ityp,i)
             enddo
           endif
           if(ityp.eq.0) then 
                  write(6,*) 'bigdia  appel qhef',ityp,i,ncf,netat
               i=1
             call extr(ncf,netat,cv,hcv,eold,qictot,ityp,i) 
           endif
         else if(qwb)then
           write(6,*)' vecteurs propres matrice davidson'
           netat=metat2
           call impression_vecteurs(ncf,netat,cv)
           endif
        endif   
        stop
      endif
      
c     Start Davidson diagonalisation

      ema=-100.
c        call prt(h,1,ncper)
        nca=ncper
          netat=0
      call model(h)
         if(qwa)call final
      write(6,989)
989   format(1x,13('*'),/,' iterations ',/,
     *13('*'),/)
          write(6,*)'starting matrix',metat2,' for',netat,' roots'
         if(netat.lt.ncper)then
         ewm=0.d0
         do190i=netat+1,ncper
190       ewm=ewm+ew(i)
          ewm=ewm/float(ncper-netat)
         else
         ewm=ew(metat2)+1.
         endif
                  nca=metat2
c         call prt(h,1,nca)
                  call model(h)
c          write(6,6005)(ew(i),i=1,nca)
c          call final
       open(unit=10,status='unknown',access='direct',recl=nvsi*8,
     *    file='BD_SCRATCH1')
         open(unit=11,status='unknown',access='direct',recl=nvsi*8,
     *    file='BD_SCRATCH2')
         write(6,*) 'bigdia appel getsec metat2',metat2
        call getsec(hcv,cv)
        do830i=1,metat2
        cw(i)=0.d0
        cw(i+netat)=0.d0
        t=0.d0
           tt=0.d0
        do820k=1,ncper
        tu=cser(k,i)*cser(k,i)
        t=t+tu
820     tt=tt+tu*cv(k)
        ew(i)=tt/t
                   cw(i)=1.d0-sqrt(t)
830     continue
c                 write(6,*)(cw(i),i=1,metat2)
c                 write(6,*)(ew(i),i=1,netat)
         do195 i=1,ncper
          qz(i)=.false.
195       cv(i)=ewm
           irhf=0
          do198i=1,ncf,nvsi
         ii=i+nvsi-1
          irhf=irhf+1
198       write(10,rec=irhf)(cv(j),j=i,ii)
             call trf(hcv,10,irhf)
      do402jit=1,njter
       nca=metat2
       mca=nca
       ircf=0
c
c   iterations davidson
c
      do400 nit=1,niter
c       do253j=1,mca
c253    write(6,254)(cold(j,i),i=1,netat)
254     format(1x,4f12.4)
        write(6,*)' ** Start Davidson : iteration n° ',nit
        write(6,*) ' ** ',mca,netat,metat2
        do230i=mca+1,mca+netat
230      qz(i-metat2)=qz(i-metat2-netat)
       nca=nca+netat
c      write(6,*)nit,(qz(i),i=1,nca-metat2)
       if(nca.gt.ncpm)then
        write(6,*)' la dimension de h effectif excede les capacites'
         nca=nca-netat
         goto401
       endif
         ij0=(mca+1)*mca/2+1
         ij=nca*(nca+1)/2
       do240i=ij0,ij
240       h(i)=0.d0
         m=nvsi
       if(mca.gt.metat2)then
       do255j=1,ncf
       do251ie=1,netat
        tc=0.d0
        t=0.d0
        ji=j-ncf
       do250 i=mca-netat+1,mca
       ji=ji+ncf
       t=t+cold(i,ie)*cv(ji)
250    tc=tc+cold(i,ie)*hcv(ji)
251     work(ie)=t*eold(ie)-tc
        ji=j-ncf
        do252i=1,netat
        ji=ji+ncf
252     cv(ji)=work(i)
255     continue
        j2=nvsi
         irh=1+(ncf-1)/nvsi
         irc=0
        do530i=1,netat
        j0=1
275     j1=j2+1
        if(j1.gt.nvsi)then
        j1=1
        irh=irh+1
        read(10,rec=irh)va
        endif
        j2=min(ncf-j0+j1,nvsi)
        do285ie=1,netat
         fac=-cold(i,ie)
               jj=ncf*ie-ncf+j0-j1
        do280j=j1,j2
280     cv(jj+j)=cv(jj+j)+fac*va(j)
285     continue
        j0=j0+j2-j1+1
        if(j0.lt.ncf)goto275
530     continue
        j2=nvsi
        do520i=metat2+1,mca-netat
        j0=1
260     j1=j2+1
        if(j1.gt.nvsi)then
         j1=1
        irh=irh+1
        irc=irc+1
        read(10,rec=irh)va
        read(11,rec=irc)vb
        endif
        j2=min(ncf-j0+j1,nvsi)
        do270ie=1,netat
         fac=-cold(i,ie)
          jj=ncf*ie-ncf+j0-j1
         fa=-fac*eold(ie)
        do265j=j1,j2
265     cv(jj+j)=cv(jj+j)+fac*va(j)+fa*vb(j)
270     continue
          j0=j0+j2-j1+1
        if(j0.lt.ncf)goto260
        if(mod(i-metat2,netat).eq.0)j2=nvsi
520     continue
        elseif(qrvec)then
        read(58) ncfx,netx,(cv(i),i=1,ncf*netat)
        else if(jit.eq.1)then
        do310j=1,netat*ncf
310     cv(j)=-hcv(j)
        endif
                 do295j=1,ncper
            ji=j
                 do294ie=1,netat
                 t=0.d0
                 do290i=1,netat
290      t=t+cold(i,ie)*cser(j,i)
       cv(ji)=cv(ji)+t*eold(ie)
        ji=ji+ncf
294      continue
295      continue
        j2=nvsi
        irh=0
         j0=1
315     j1=j2+1
        if(j1.gt.nvsi)then
        irh=irh+1
                  j1=1
        read(10,rec=irh)va
        endif
        j2=min(ncf-j0+j1,nvsi)
         if(.not.(qrvec.or.jit.gt.1).or.mca.gt.metat2)then
        do325ie=1,netat
          ezer=-eold(ie)
         jj=ncf*ie-ncf+j0-j1
c          if(qz(mca-metat2+ie))then
c          do321j=j1,j2
c321        cv(jj+j)=0.d0
c          else
          t=0.d0
        do320j=j1,j2
320     cv(jj+j)=cv(jj+j)/(ezer+va(j))
c        endif
325     continue
        endif
        if(mca.le.metat2)then
        do855ie=1,netat
         jj=ncf*ie-ncf+j0-j1
        t=0.d0
        tt=0.d0
        ja=j1
         if(j0.le.ncper)ja=ncper+1
         do850j=ja,j2
        t=t+hcv(jj+j)*hcv(jj+j)/(ew(ie)-va(j))
        tt=tt+hcv(jj+j)*hcv(jj+j)/(eold(ie)-va(j))
850     continue
        cw(ie)=cw(ie)+t
        cw(ie+netat)=cw(ie+netat)+tt
855     continue
        endif
                   j0=j0+j2-j1+1
        if(j0.lt.ncf)goto315
         if(mca.eq.metat2)then
         write(6,*)
         write(6,*)' second ordre '
         write(6,*)' d bary ',(ew(i),i=1,netat)
         write(6,*)' e bary ',(cw(i),i=1,netat)
         write(6,*)' e v p ',(cw(i+netat),i=1,netat)
         write(6,*)
         endif
6006  format(1x,14f9.4)
        call mproj(cv)
       call oproj(cv,hcv,mca-metat2)
       call low(cv,ncf)
      call hmcv(cv,hcv)
        call trf(hcv,10,irhf)
        call trf(cv,11,ircf)
            ik=mca*(mca+1)/2
            ji=-ncf
          do365i=1,netat
            ji=ji+ncf
          do355k=1,metat2
             t=0.d0
c       if(qz(i+mca-metat2))goto351
          do350j=1,ncper
350        t=t+hcv(ji+j)*cser(j,k)
351           h(ik+k)=t
355        continue
              ik=ik+mca
         jk=0
          do365k=1,i
          ik=ik+1
             t=0.d0
c       if(qz(k+mca-metat2).or.qz(i+mca-metat2))then
c        if(k.eq.i)t=100.d0
c        jk=jk+ncf
c        goto361
c       endif
          do360j=1,ncf
           jk=jk+1
360       t=t+cv(ji+j)*hcv(jk)
361           h(ik)=t
365        continue
           irc=0
           j2=nvsi
          iie=mca*(mca+1)/2
        do540i=metat2+1,mca
        j0=1
375     j1=j2+1
        if(j1.gt.nvsi)then
        j1=1
        irc=irc+1
        read(11,rec=irc)vb
        endif
        j2=min(ncf-j0+j1,nvsi)
          iie=mca*(mca+1)/2+i
        do385ie=1,netat
         jj=ncf*ie+j0-ncf-j1
          t=0.d0
c          if(qz(i-metat2).or.qz(ie+mca-metat2))goto381
        do380j=j1,j2
380      t=t+hcv(jj+j)*vb(j)
381         h(iie)=h(iie)+t
         iie=iie+ie+mca
385     continue
          j0=j0+j2-j1+1
        if(j0.lt.ncf)goto375
        if(mod(i-metat2,netat).eq.0)j2=nvsi
540     continue
        if(nit.eq.1)call prt(h,1,nca)
      call model(h)
       write(6,*)(eold(iwi),iwi=1,netat)
      it=0
      do110 k=1,netat
        iperf(k)=abs(cold(k+mca,k))*float(ites(k))
c        qz(k+mca-metat2)=iperf(k).lt.1
      if(iperf(k).gt.it)then
       it=iperf(k)
      endif
110       continue
          write(6,*)' iter',nit,'  perf',(iperf(k),k=1,netat)
      if(it.lt.1.and.nit.gt.1)goto6001
          write(6,6004)(eold(k),k=1,netat)
        mca=mca+netat
400        continue
401         call wvec(cv,hcv,nca-metat2,.true.)
402         continue
         write(6,*)' non converge voir niter'
6001       write(6,6005)(ew(i),i=1,nca)
6005   format(' final result, energies',/,(1x,12f10.4))
6003   write(6,6004)(eold(i),i=1,netat)
6004   format(' ref space energies',(1x,8f14.6))
      typener='energies var.'//title
c
c     stockage des energies
c
      call stkener(typener,eold,netat)
c
c  ecriture des vecteurs 
c
        if(qhef.or.qwvec.or.qws)then
           call wvec(cv,hcv,nca-metat2,.false.)
           if(qwvec)then
            ideb=1
            ifin=ncf*netat
            ncf59=ncf
c     recombinaison sur les determinants
            if (qdcple) then
             rac2p=1.d0/dsqrt(2.0d0)
             read(4)
             read(4) isymat,(qef(i),i=1,ncf)
             if (isymat.eq.1) rac2pm=rac2p
             if (isymat.eq.3) rac2pm=-rac2p
             do i=1,ncf
                write(6,*) qef(i)
                 enddo
             netncf=netat*ncf*2
             inda=0
             indb=netncf/2
             do j=netat,1,-1
               do i=j*ncf,((j-1)*ncf+1),-1
                    indc=ncf
                if (.not.qef(indc)) then
                   cv(netncf-inda)=cv(indb)
                   inda=inda+1
                    else
                   cv(netncf-inda)=cv(indb)*rac2p
                   inda=inda+1
                   cv(netncf-inda)=cv(indb)*rac2pm
                   inda=inda+1
                    endif
                    indc=indc-1
                indb=indb-1
             enddo
             enddo
             ifin=netncf
             ideb=netncf-inda+1
             ncf59=inda/netat
            endif
             write(6,*) 'ifin,ideb,ncf,ncf*netat,ncf59',ifin,
     &                       ideb,ncf,ncf*netat,ncf59
c
c      ecriture des vecteurs et valeurs propres completes sur l unite 59
c
              rewind 59
              write(59) ncf59,netat,(cv(i),i=ideb,ifin),
     *            (eold(i)+escf,i=1,netat)
           endif
c
c     Calcul de l energie totale et des corrections
c     Hamiltonien effectif
c
           if(qhef) then
           if(ityp.ne.0) then 
             do i =1,3
                  write(6,*) 'bigdia  appel qhef',ityp,i,ncf,netat
              call extr(ncf,netat,cv,hcv,eold,qictot,ityp,i)
             enddo
           endif
           if(ityp.eq.0) then 
                  write(6,*) 'bigdia  appel qhef',ityp,i,ncf,netat
               i=1
             call extr(ncf,netat,cv,hcv,eold,qictot,ityp,i) 
           endif
         else if(qwb)then
           write(6,*)' vecteurs propres matrice davidson'
           netat=metat2
           call impression_vecteurs(ncf,netat,cv)
           endif
       endif
       stop
      end
