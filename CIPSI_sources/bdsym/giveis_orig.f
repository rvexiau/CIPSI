      subroutine giveis(n,nv,nd,a,b,root,vect,ierr)
      implicit real*8(a-h,o,p,r-z),logical*4(q)
      dimension a(1),vect(1)
      common/info/ncf,ncper,nca,qnorm,qwa,qwb,mbuf,metat2,netat,niter
     *,qwvec,qhef,ntrsy,qrvec,seuil,qws,qsuis,qsuip,qener,njter,qgiv
     *,qictot,qsui04,qdcple
      if (qgiv) then
         nev=n
         nvec=n
         ijsup=0
         do 5  i=1,n
         do 5  j=i,n
      ijsup=ijsup+1
      ijinf=j*(j-1)/2+i
5     vect(ijsup)=a(ijinf)
      nij=n*(n+1)/2
      do 8 ijsup=1,nij
8     a(ijsup)=vect(ijsup)
      call givens(a,root,vect,n,nev,nvec)
      else
         call jacscf(a,vect,root,n,-1,1.d-12)
      endif
      return
      end
