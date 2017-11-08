      subroutine tmono(ipass,nveci,nvecf,ch,csc,yfi50)
      implicit real*8 (a-h,o-x,z),logical*1 (y)
      include 'pshf.prm'
      dimension csc(2*doa*doa)
      dimension ch(doa*(doa+1)/2)
      common norb,nsym,its(nsymz,nsymz),itsym(doa),itsyv(2*doa),
     1 nvec,isymp(nsymz),mijkl,nrec,mblk
      common isydeg(nsymz,nsymz),ntrsy,istop,iarwr(2*doa),
     1 nocsy(nsymz),ymono,yprt,yprtv,yatom,ywvect,ynocsy
      common /om/ c(doa,2*doa),nbp(nsymz),nbo(nsymz),ind(nsymz+1),
     y ivecsy(2*doa*nsymz)
      common /hondo/ nt,ishell(doa),kloc(doa),ktyp(doa),
     1 numc(doa),newsh(doa,nsymz),iptr(3,nsymz),idtr(6,nsymz),
     2 iftr(10,nsymz),igtr(15,nsymz),yhondo,yprth
c
      nij=norb*(norb+1)/2
      nnorb=norb*norb
c lecture des integrales monoelectroniques
      rewind 2
      read(2)
      read(2)
      read(2) (ch(ij),ij=1,nij)
c
c
c     relecture des coeff. sur la file 10 ,si fi50=t
c
      if(yfi50) then
      if(ipass.eq.1)then
      rewind 10
      call readb (csc,nnorb,10)
      else
      call readb(csc(nnorb+1),nnorb,10)
      end if
      end if
c
c     write(6,*) 'fi50',yfi50
c     write(6,*) 'ch',(ch(i),i=1,nij)
c     write(6,*) 'csc',(csc(i),i=1,nnorb)
      do 200 i=nveci,nvecf
      ni=(i-1)*norb
      do 200 j=1,norb
      c(j,i)=0.d0
      do 200 k=1,norb
      if(k.ge.j) kj=k*(k-1)/2+j
      if(k.lt.j) kj=j*(j-1)/2+k
  200 c(j,i)=c(j,i)+csc(ni+k)*ch(kj)
      ij=0
      do 300 i=nveci,nvecf
      do 300 j=nveci,i
      nj=(j-1)*norb
      ij=ij+1
      ch(ij)=0.d0
      do 300 k=1,norb
  300 ch(ij)=ch(ij)+c(k,i)*csc(nj+k)
c
      if(yfi50)write(50) (ch(i),i=1,ij)
          if(.not.yprt) return
      write(6,9998)
      jf=0
      do 910 i=1,norb
      ji=jf+1
      jf=jf+i
  910 write(6,9997) (ch(j),j=ji,jf)
 9998 format(/,' integrales monoelectroniques (t+v) dans la base molecul
     &aire')
 9997 format(12f15.6)
      return
      end
