      subroutine  trf(hcv,ifi,irhf)
      implicitreal*8(a-h,o-p,r-z),logical*4(q)
      integer smax      
      include 'bd.prm'
      common/gest/spert,eold(netm),cold(ncpm,netm),ind(ncpm),qz(ncpm)
      common/info/ncf,ncper,nca,qnorm,qwa,qwb,mbuf,metat2,netat,niter
     *,qwvec,qhef,ntrsy,qrvec,seuil,qws,qsuis,qsuip,qener,njter,qgiv
     *,qictot,qsui04,qdcple
      common/serre/cser(nczm,nczm),va(nvsi),vb(nvsi),ites(netm)
      dimension hcv(nvsi)
          do198i=1,ncf*netat,nvsi
        ii=i+nvsi-1
          irhf=irhf+1
198       write(ifi,rec=irhf)(hcv(j),j=i,ii)
            return
            end
