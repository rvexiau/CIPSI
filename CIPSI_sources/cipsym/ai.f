      function ai(i,j,k,l)                                       
      include 'pshf.prm'
      integer*2 ij0 
      real*4 bijkl                                      
      real*8 ai
      common/ij1/ij0(doa,doa),ntttt(nsymz)                           
      common/bij/bijkl(kget)
      if(j.ge.i) then
         ij=ij0(i,j)
         if(i.eq.j) then
            isym=1
            else
            isym=ij0(j,i)
         endif
         else
         ij=ij0(j,i)
         isym=ij0(i,j)
      endif
      if(l.ge.k) then
         kl=ij0(k,l)
         if(k.eq.l) then
            jsym=1
            else
            jsym=ij0(l,k)
         endif
         else
         kl=ij0(l,k)
         jsym=ij0(k,l)
      endif
      if(isym.ne.jsym) then
         ai=0.d0
         return
         else
         if(ij.le.kl) then
            ijkl=ij+(kl*kl-kl)/2+ntttt(isym)-1
            else
            ijkl=kl+(ij*ij-ij)/2+ntttt(isym)-1
         endif
      endif       
      ai=bijkl(ijkl)
      return
      end
