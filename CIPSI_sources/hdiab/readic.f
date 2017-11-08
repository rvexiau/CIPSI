      subroutine readic(v,h,ne,nd,tr,pr,nv,ncf,irvec,irdet,ic)
      implicit real*8(a-h,o,p,r-z),logical*1(q)
      include 'pshf.prm' 
      integer*4 nttt(2*nexz)
      integer*4 ne,nd,tr,pr
      dimension v(metz*ndetz),ne(ndetz),nd(ndetz),tr(ntpz),pr(ntpz)
      dimension h(metz,metz)
c
c   Lecture vecteurs d ic et determinants
c   ic=1   lecture moyen+bd
c   ic=0   lecture CIPSI uniquement
c
      do i=1,metz
      do j=1,metz
        h(i,j)=0.d0
      enddo
      enddo
      if(ic.eq.1) then
      write(6,*) 'lecture dets moyen et vecteurs bd'
      write(6,*) 'files:  vec ',irvec, ' det ',irdet
c
c    lecture moyen-davidson
c
       read(irvec) ncf,nv,((v(i+(j-1)*ncf),i=1,ncf),j=1,nv),
     *              (h(i,i),i=1,nv)
       write(6,*) 'ncf ',ncf,' nvec ',nv
       ndf=0
       do idet=1,ncf
       nd(idet)=ndf
       read(irdet,end=6341)nec
       ne(idet)=nec
       if(nec.ne.0)read(irdet,end=6341)(nttt(j),j=1,nec+nec)
           do j=1,nec
                ndf=ndf+1
                tr(ndf)=nttt(j)
                pr(ndf)=nttt(j+nec)
           enddo
       enddo
       
      write(6,*) ' determinants lus ',ncf
      else
c
c    lecture cipsi
c
      read(irvec) norb,noca,nocb,nv,ncf,isz,qion,ii,ik,
     * (ne(j),nd(j),j=1,ncf),(tr(j),pr(j),j=1,ii),
     * (v(j),j=1,ik)
       endif
c  
c  ecriture vecteurs
c
      write(6,*) 'vecteurs'
      call wrtvec(v,nv,ncf)
       write(6,*) 'fin lecture IC'
       return
6341   write(6,*) ' erreur lecture sur file moyen ', n1fim
       stop
       end
