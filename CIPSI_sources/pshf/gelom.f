      subroutine gelom(v,e,da,db,fa,fb,ia,in,uhf)
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      logical uhf
c
      common/pseudo/qa(doa),qb(doa),ngla(doa),nglb(doa),numgla,numglb
      common/iofile/ir,iw,ip,is,iq,ih,iv
      common/infoa/nat,ich,mul,numscf,nx,ne,na,nb
      dimension v(*),e(*),fa(*),fb(*),da(*),db(*),ia(*),in(*)
      dimension g(doa)
 9999 format(//' interdit de geler des om couches ouvertes..',i3,f5.2,
     *'    pas de gel')
 9998 format(//,' gel des om. construction de h+sig(j-k) et reclassement
     *des vecteurs')
 9997 format(//' vecteurs apres reclassement')
 9996 format(/'  energie des om gelees',f20.8)
      if(numgla.eq.0) return
      write(iw,9998)
      nxsq=numscf*numscf
      call reada(v,nxsq,7)
      do 100 l=1,numgla
      k=ngla(l)
      if(qa(k).eq.0.) go to 100
      if(qa(k).eq.2.) go to 80
      write(iw,9999) k,qa(k)
      return
   80 nk=in(k)
      ij=0
      do 90 i=1,numscf
      do 90 j=1,i
      ij=ij+1
   90 fa(ij)=fa(ij)+v(nk+i)*v(nk+j)*2.d0
  100 continue
      rewind is
      call hstar(da,fa,db,fb,ia,uhf)
      call symh(fa,da,ia)
      call reada(fa,nx,3)
      call addmat(fa,da,da,nx)
      call wrtda(da,nx,3)
      call addmat(fa,da,da,nx)
      call reada(v,nxsq,7)
      dum=0.d0
      do 140 i=1,numgla
      k=ngla(i)
      if(qa(k).ne.2.) go to 140
      nk=in(k)
      do 150 j=1,numscf
      do 150 l=1,j
      jl=ia(j)+l
      fx=da(jl)
      if(j.ne.l) fx=fx*2.0d0
      dum=dum+fx* v(nk+j)*v(nk+l)
  150 continue
  140 continue
      write(iw,9996) dum
      ij=0
      do 160 i=1,numgla
      k=ngla(i)
      g(i)=e(i)
      do 160 j=1,numscf
      ij=ij+1
  160 da(ij)=v(in(k)+j)
      do 200 i=1,numgla
      k=ngla(i)
c
c
      ji=k-i+1
      jf=numscf-i
      if(ji.gt.jf) go to 200
      do 210 j=ji,jf
      e(j)=e(j+1)
      do 210 l=1,numscf
  210 v(in(j)+l)=v(in(j+1)+l)
      j=numscf-numgel+i
c
c
  200 continue
      ik=0
      ij=numscf*(numscf-numgla)
      do 215 i=1,numgla
      e(numscf-numgla+i)=g(i)
      do 215 j=1,numscf
      ij=ij+1
      ik=ik+1
  215 v(ij)=da(ik)
      write(iw,9997)
      call vecout(v,e,in,numscf,numscf)
      call wrtda(v,nxsq,7)
      return
      end
