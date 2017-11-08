!
!     DB (06/2014) : 
!       Appel à 'giveis' remplacés par un appel à 'diagonaliser'
!

      subroutine  low(cv,ncf)
      implicitreal*8(a-h,o-p,r-z),logical*4(q)
      include 'bd.prm'
      common/gest/spert,eold(netm),cold(ncpm,netm),ind(ncpm),qz(ncpm)
      common/info/mcf,ncper,nca,qnorm,qwa,qwb,mbuf,metat2,netat,niter
     *,qwvec,qhef,ntrsy,qrvec,seuil,qws,qsuis,qsuip,qener,njter,qgiv
     *,qictot,qsui04,qdcple
      common/serre/cser(nczm,nczm),va(nvsi),vb(nvsi),ites(netm)
      dimension c(ncpm*9),cv(ncf*netat),qeli(netm)
       common/wry/ work(ncpm*(ncpm+1)/2),cw(ncpm*ncpm),ew(ncpm)
      equivalence(c(1),va(1))
      if(netat.eq.1)return
       liter=0
1     continue
       liter=liter+1
       if(liter.gt.3)then
       write(6,*)' precision douteuse low  apres 2iter'
       return
       endif
         kj=-ncf
               kl=0
         do250k=1,netat
        qeli(k)=.false.
         lj=0
         kj=kj+ncf
         do250l=1,k
         t=0.d0
         do220j=1,ncf
         lj=lj+1
220      t=t+cv(kj+j)*cv(lj)
         kl=kl+1
         work(kl)=t
         if(work(kl).lt.0.9d0)qeli(k)=.true.
250      continue
       nd=netat
       !call giveis(nd,nd,nd,work,c,ew,cw,ierr)
       call diagonaliser( nd, work, ew, cw,nd )
       neli=0
      t=0.d0
      do 130 i=1,nd
      t=t+ew(i)
       if(ew(i).gt.1.d-8)then
      ew(i)=1.d0/dsqrt(ew(i))
          else
         neli=neli+1
         ew(i)=0.d0
         endif
130      continue
c      write(6,*)t,neli
      do 150 i=1,nd
      do150j=1,i
       t=0.d0
       jk=j-nd
       ik=i-nd
      do 140 k=1,nd
      jk=jk+nd
      ik=ik+nd
140     t=t+cw(ik)*ew(k)*cw(jk)
       c((i-1)*nd+j)=t
       c((j-1)*nd+i)=t
150   continue
      do 180 ip=1,ncf
      ik=0
      do 160 i=1,nd
      t=0.d0
       ki=ip-ncf
      do 155 k=1,nd
       ki=ki+ncf
         ik=ik+1
155   t=t+cv(ki)*c(ik)
160   ew(i)=t
          ki=ip-ncf
          kk=0
      do 170 i=1,nd
       kk=kk+1
        ki=ki+ncf
        cv(ki)=ew(kk)
170     continue
180      continue
         ii=0
        do195i=1,nd
        t=0.d0
        if(qeli(i))goto195
        do190k=1,nd
        ii=ii+1
190     t=t+c(ii)
        t=t/c(ii-nd+i)
        if(abs(t).lt.1.d-5)then
        write(6,*) ' pb suppose low prec',t
        goto1
         endif
195     continue
         return
         end
