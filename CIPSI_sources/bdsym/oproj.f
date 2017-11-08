      subroutine oproj(c,cc,mm)
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
      do11k=1,netat
      rec=0.d0
      do10j=1,ncf
      jk=jk+1
10    rec=rec+c(jk)*c(jk)
11    ov(k)=rec
c        write(6,61)(c(j),j=1,ncf*netat)
c61      format(   8f10.4)
c      write(6,*)' norme',(ov(k),k=1,netat)
      do100iter=1,25
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
c         if(qz(k))goto130
        ji=0
      do30i=1,netat
      t=0.d0
      do20j=1,ncf
20    t=t+c(ji+j)*cc(j)
      ov(i)=ov(i)-t*t
      do 25 j=1,ncf
       ji=ji+1
25    c(ji)=c(ji)-t*cc(j)
30    continue
130         continue
        jk=-ncf
c                 write(6,*)' norm2',(ov(k),k=1,netat)
       qit=.false.
        kj=-ncf
       do90k=1,netat
       kj=kj+ncf
       rec=ov(k)
      if(rec.gt.1.d-20) then
      rec=1.d0/sqrt(rec)
        else
c      elseif(.not.qz(mm+k))then
      write(6,*)' elimine',k,rec
c        qz(mm+k)=.true.
      rec=0.d0
      endif
      tt=0.d0
      do60j=1,ncf
      t=c(kj+j)*rec
      tt=t*t+tt
60     c(kj+j)=t
      ov(k)=tt
      if(rec.eq.0.d0)goto90
      if(abs(tt-1.d0).gt.1.d-12)then
      qit=.true.
       rec=1.d0/dsqrt(tt)
        endif
90    continue
cc       write(6,61)(ov(jj),jj=1,netat)
c       write(6,61)(c(jj),jj=1,ncf*netat)
      if(.not.qit)return
       write(6,*)' oproj***   ',(ov(k),k=1,netat)
100    continue
      write(6,*)' vecteur ',k,' norme',tt,'malgre'
      return
      end
