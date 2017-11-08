      subroutine wvec(c,cc,mm,qjit)
      implicitreal*8(a-h,o-p,r-z),logical*4(q)
      include 'bd.prm'
      common/gest/spert,eold(netm),cold(ncpm,netm),ind(ncpm),qz(ncpm)
      common/info/ncf,ncper,nca,qnorm,qwa,qwb,mbuf,metat2,netat,niter
     *,qwvec,qhef,ntrsy,qrvec,seuil,qws,qsuis,qsuip,qener,njter,qgiv
     *,qictot,qsui04,qdcple
      common/serre/cser(nczm,nczm),va(nvsi),vb(nvsi),ites(netm)
       common/wry/ work(ncpm*(ncpm+1)/2),cw(ncpm*ncpm),ew(ncpm)
      dimension c(ncf*netat),cc(ncf),ov(ncpm)
       jk=0
       do100i=1,netat
       do100j=1,ncf
       jk=jk+1
100    c(jk)=0.d0
       j0=0
        irc=0
        j2=nvsi
      do130k=1,mm
110     j1=j2+1
        if(j1.gt.nvsi)then
        j1=1
         irc=irc+1
        read(11,rec=irc)va
        endif
        j2=min(ncf-j0+j1-1,nvsi)
        do120j=j1,j2
         j0=j0+1
120      cc(j0)=va(j)
          if(j0.lt.ncf)goto110
        if(mod(k,netat).eq.0)j2=nvsi
           j0=0
        ji=0
      do30i=1,netat
      t=cold(metat2+k,i)
      do 25 j=1,ncf
       ji=ji+1
25    c(ji)=c(ji)+t*cc(j)
30    continue
130         continue
        jk=-ncf
      do50k=1,netat
       jk=jk+ncf
      do50i=1,metat2
      t=cold(i,k)
      do 45 j=1,ncper
45    c(jk+j)=c(jk+j)+t*cser(j,i)
50    continue
       if(qjit)return
          if(qwb.or.qws)then
          ki=0
          kj=-ncf
         do60k=1,netat
          kj=kj+ncf
         do60i=1,ncper
          ki=ki+1
60       cw(ki)=c(kj+i)
          nca=ncper
         metat2=netat
      write(6,987)
987   format(/,18('*'),/,' resultats finaux ',/,18('*'),/)
      write(6,*) ' ***  norme des etats ***'
      write(6,*) ' etat       norme dans S      norme total'       
      do mm=1,netat
        dnorme=0d0
        dnormes=0d0
        do i=1,nca
          dnorme= dnorme + c((mm-1)*ncf+i)*c((mm-1)*ncf+i)
        enddo       
        dnormes=dnorme
        do i=nca+1,ncf
          dnorme= dnorme + c((mm-1)*ncf+i)*c((mm-1)*ncf+i)
        enddo   
        write(6,'(1x,i3,4x,f15.8,2x,f15.8)') mm,dnormes,dnorme
      enddo  
      write(6,*)
      write(6,*)' composantes des vecteurs propres dans s'
             if(qwb)call final
      if(qws)write(55)ncper,netat,(cw(k),k=1,ki),(eold(k),k=1,netat)
           endif
      return
      end
