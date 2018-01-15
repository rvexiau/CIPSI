      subroutine storic(iwant,io)
      implicit real*8(a-h,o-p,r-z),logical*4(q)
      integer smax      
      include 'bd.prm'
      common/gest/spert,eold(netm),cold(ncpm,netm),ind(ncpm),qz(ncpm)
      common/info/ncf,ncper,nca,qnorm,qwa,qwb,mbuf,metat2,netat,niter
     *,qwvec,qhef,ntrsy,qrvec,seuil,qws,qsuis,qsuip,qener,njter,qgiv
     *,qictot,qsui04,qdcple
      common/serre/cser(nczm,nczm),va(nvsi),vb(nvsi),ites(netm)
      common/dyn/d(kget)
      
c          iget=ia01as(d,iwant,8)
c           if(iget.lt.0)then
            if(iwant.gt.kget)then
              write(6,*)' taille demandee',iwant,'  superieure a la
     *taille prevue(kget)',kget
              stop
            endif
         iget=1
      write(6,*) 'storic'
      write(6,*) 'qnorm',qnorm
c        if(qnorm)call norma(d(iget))
       mpos=max(ncf*netat,ncper*ncper)
      write(6,*) 'storic'
      write(6,*) 'ncf,netat,ncper,iget'

      write(6,*) ncf,netat,ncper,iget
      call bigdia(d(iget),d(mpos+iget))
      write(6,*) 'end of bigdia' 
      
      return
      end
