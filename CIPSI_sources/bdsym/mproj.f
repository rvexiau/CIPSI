      subroutine mproj(c)
      implicitreal*8(a-h,o-p,r-z),logical*4(q)
      integer smax      
      include 'bd.prm'
      common/gest/spert,eold(netm),cold(ncpm,netm),ind(ncpm),qz(ncpm)
      common/info/ncf,ncper,nca,qnorm,qwa,qwb,mbuf,metat2,netat,niter
     *,qwvec,qhef,ntrsy,qrvec,seuil,qws,qsuis,qsuip,qener,njter,qgiv
     *,qictot,qsui04,qdcple
      common/serre/cser(nczm,nczm),va(nvsi),vb(nvsi),ites(netm)
       common/wry/ work(ncpm*(ncpm+1)/2),cw(ncpm*ncpm),ew(ncpm)
      dimension c(ncf*netat),ov(ncpm)
       jk=0
      do11k=1,netat
      rec=0.d0
      do10j=1,ncper
10    rec=rec+c(jk+j)*c(jk+j)
        jk=jk+ncf
11    ov(k)=rec
                  qit=.false.
             if(rec.gt.0)then
             tre=1.d-10/sqrt(rec)
             else
             tre=1.d-10
             endif
      do100iter=1,25
       qit=.false.
        jk=-ncf
      do50k=1,netat
       jk=jk+ncf
      do50i=1,metat2
      t=0.d0
      do40j=1,ncper
40    t=t+cser(j,i)*c(jk+j)
      ov(k)=ov(k)-t*t
      do 45 j=1,ncper
45    c(jk+j)=c(jk+j)-t*cser(j,i)
       if(abs(t).gt.tre)qit=.true.
50    continue
      if(.not.qit)return
100    continue
       write(6,*)' mproj***   ',(ov(k),k=1,netat)
      return
      end
