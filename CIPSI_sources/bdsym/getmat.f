      subroutine  getmat(h,cv,hcv,io)
      implicitreal*8(a-h,o-p,r-z),logical*4(q)
      include 'bd.prm'
      common/gest/spert,eold(netm),cold(ncpm,netm),ind(ncpm),qz(ncpm)
      common/info/ncf,ncper,nca,qnorm,qwa,qwb,mbuf,metat2,netat,niter
     *,qwvec,qhef,ntrsy,qrvec,seuil,qws,qsuis,qsuip,qener,njter,qgiv
     *,qictot,qsui04,qdcple
      common/serre/cser(nczm,nczm),va(nvsi),vb(nvsi),ites(netm)
      dimension h(1),cv(1),tab(lsize)
      dimension d(lsize),c(ncf*netm),b(ncf*netm),hcv(1)
      integer*4 itab(lsize),jtab(lsize)
      save lpos,itab,jtab,d,last
c
c    lecture de la partie  la matrice _mat de moyen qui
c    sert a la diagonalisation initiale ( dimension ncper)
c
c    h contient la demi matrice inferieure
c    hcv contient la matrice entiere
c
       ij=0
       do10i=1,ncper
       do10j=1,i
       hcv(ncper*(j-1)+i)=0.d0
       hcv(ncper*(i-1)+j)=0.d0
      ij=ij+1
10     h(ij)=0.d0
       rewind(io)
11       read(io,end=210)itab,jtab,d
       do iii=1,lsize
       enddo
       if(itab(lsize).lt.0)then
        last=jtab(lsize)
       else
         last=lsize
       endif
       do17l=1,last
        i=itab(l)
        j=jtab(l)
        if(i.gt.ncper)goto25
      h(i*(i-1)/2+j)=d(l)
      hcv(i+ncper*(j-1))=d(l)
      hcv(j+ncper*(i-1))=d(l)
17    continue
            lpos=l
          if(itab(lsize).lt.0)then
          if(ncper.eq.i)goto25
          return
           endif
          goto11
25        continue
          lpos=l
                ij=0
      do 15 i=1,ncper
      ij=ij+i
            h(ij)=h(ij)+h(ij)
            hcv(i+(i-1)*ncper)=h(ij)
15          continue
      ii=0
         return
c ***********************************************************
          entry getsec(cv,hcv,io)
          write(6,*) ' debut getsec ncf,ncper,netat,lpos,last'
          write(6,*) ncf,netat,ncper,lpos,last
c
c   lecture de diagonale de la matrice pour calcul de la perturbation
c   hcv contient la diagonale
c
            do28j=1,ncf*netat
28          cv(j)=0.d0
                   kj=-ncf
                   do29k=1,netat
          ij=0
          kj=kj+ncf
         do27j=1,ncper
         t=0.d0
         do26i=1,ncper
          ij=(j-1)*ncper+i
26       t=t+cser(i,k)*hcv(ij)
27       cv(kj+j)=t
29       continue
            ki=-ncper
          do340k=1,ncper
          ki=ki+ncper
340       hcv(k)=hcv(ki+k)
                l1=lpos
30      do100l=l1,last
         i=itab(l)
         j=jtab(l)
         if(j.le.ncper)then
            ki=i-ncf
           do40k=1,netat
           ki=ki+ncf
           cv(ki)=cv(ki)+cser(j,k)*d(l)
40         continue
c ne marchera pas si matrice random
           elseif(j.eq.i)then
           hcv(i)=d(l)+d(l)
           endif
100    continue
          if(itab(lsize).gt.0)then
           l1=1
           read(io,end=210) itab,jtab,d
           if(itab(lsize).lt.0)then
             last=jtab(lsize)
           else
             last=lsize
           endif
             goto30
            endif
c         write(6,61)(cv(j),j=1,ncf*netat)
c61       format(8f10.4)
        return
c ***********************************************************
        entry hmcv(b,c,io)
c
c   lecture matrice
c
       rewind(io)
         k0=nca-netat-metat2
       do150i=1,ncf*netat
150    c(i)=0.d0
160       read(io,end=210)itab,jtab,d
       if(itab(lsize).lt.0)then
        last=jtab(lsize)
       else
         last=lsize
       endif
       do200l=1,last
        i=itab(l)
        j=jtab(l)
         kk=-ncf
         do190k=1,netat
          kk=kk+ncf
c          if(qz(k+k0))goto190
       c(kk+i)=c(kk+i)+d(l)*b(kk+j)
        c(kk+j)=c(kk+j)+d(l)*b(kk+i)
190     continue
200    continue
                if(itab(lsize).gt.0)goto160
      return
210         write(6,*)' fin fichier sur1 i,j',i,j
c **************************************************************
      entry norma(tab,io)
         qend=.false.
      ij=0
      rewind2
      rewind(io)
                   ll=0
      kbuf=mbuf
        l=kbuf+1
      do400i=1,ncf
              do500j=1,i
              if(ll.ge.lsize)then
              write(io)itab,jtab,d
              ll=0
              endif
                l=l+1
                if(l.gt.kbuf)then
       if(mbuf.eq.0)kbuf=i
                l=1
                call rea2(tab,kbuf,qend)
                  if(qend)goto590
                endif
              if(abs(tab(l)).lt.1.d-9)goto500
              ll=ll+1
              itab(ll)=i
              jtab(ll)=j
              d(ll)=tab(l)
500           continue
             if(jtab(ll).eq.i)then
       d(ll)=.5d0*d(ll)
             else
           ll=ll+1
           itab(ll)=i
           jtab(ll)=i
           d(ll)=tab(l)
           endif
400    continue
590           if(qend)ncf=kbuf
            if(ll.eq.lsize)then
              write(io)itab,jtab,d
              ll=0
             endif
              itab(lsize)=-1
              jtab(lsize)=ll
              write(io)itab,jtab,d
       rewind(io)
       return
            end
