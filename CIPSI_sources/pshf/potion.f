      subroutine potion(f,v,e,ia,in,npi)
      implicit real*8(a-h,o-z)
      include 'pshf.prm'
      logical exch
      dimension g(doa),h(doa)
      common/restar/timlim,irest,nrec,intloc,ist,jst,kst,lst
      common/iofile/ir,iw,ip,is,iq,ih,ix
      common/infoa/nat,ich,mul,numscf,nx,ne,na,nb,zan(dc),c(3,dc)
      common/output/iop1,itol,icut,normf,normp,minp
c     common/eigen/f(doas),v(doasq),e(doa),ia(doa),in(doa)
      dimension f(*),v(*),e(*),ia(*),in(*)
 1000 format(1h ,//,10x,21(1h-),/,1h ,9x,'ionization potentials',
     #/,1h ,9x,21(1h-),/)
 1010 format(1h ,' epsilon (',i2,') = ',f6.2,' ev',4x,
     #'repolarization = ',f5.2,' ev',4x,'correlation = ',f5.2,' ev',
     #4x,'ip(l=1/2) = ',f5.2,' ev  ip(l=1) = ',f5.2,'ev')
      nxsq=numscf*numscf
      nap=na+1
      write(iw,1000)
      do 100 ii=1,npi
      ion=na+1-ii
      rewind is
      ni=nrec+4
      do 5 i=1,ni
    5 read(is)
      read(is) (f(i),i=1,nx)
      read(is) (v(i),i=1,nxsq)
      jk=0
      do 10 j=1,numscf
      do 10 k=1,j
      jk=jk+1
 10   f(jk)=f(jk)-v(in(ion)+j)*v(in(ion)+k)
      rewind is
      exch=.true.
      call hstar(exch)
      call symh
      read(is) (f(i),i=1,nx)
      ij=0
      do 20 i=1,numscf
      do 20 j=1,i
      ij=ij+1
   20 f(ij)=f(ij)+v(ij)
      read(is)
      read(is)
      read(is)
      read(is) (v(j),j=1,nxsq)
      do 60 i=1,numscf
      do 40 j=1,numscf
      dum=0.d0
      do 35 k=1,numscf
      ni=in(i)+k
      if(j.lt.k) go to 30
      jk=ia(j)+k
      go to 35
 30   jk=ia(k)+j
 35   dum=dum+f(jk)*v(ni)
 40   g(j)=dum
      do 50 j=i,numscf
      dum=0.d0
      do 45 k=1,numscf
      ni=in(j)+k
 45   dum=dum+v(ni)*g(k)
 50   v(in(i)+j)=dum
 60   continue
      do 70 i=1,numscf
      do 70 j=i,numscf
 70   f(ia(j)+i)=v(in(i)+j)
      dum=0.d0
      do 80 i=1,na
      dmult=2.d0
      if(i.eq.ion) dmult=1.d0
      do 80 j=nap,numscf
      ji=ia(j)+i
 80   dum=dum+dmult*f(ji)**2/(e(i)-e(j))
c
      eiev=27.21*e(ion)
      dumev=-27.21*dum
c
c                              -
c     correlation de la paire ii
c
      backspace is
      read(is) (v(i),i=1,nxsq)
      ij=0
      do 200 i=1,numscf
      do 200 j=1,i
      ij=ij+1
 200  f(ij)=-v(in(ion)+i)*v(in(ion)+j)
      rewind is
      exch=.false.
      call hstar(exch)
      call symh
      do 210 i=1,nx
  210 f(i)=v(i)*2.d0
      read(is)
      read(is)
      read(is)
      read(is)
      read(is) (v(i),i=1,nxsq)
      do 250 i=nap,numscf
      do 230 j=1,numscf
      dum=0.d0
      do 220 k=1,numscf
      if(j.lt.k) go to 222
      jk=ia(j)+k
      go to 220
 222  jk=ia(k)+j
 220  dum=dum+f(jk)*v(in(i)+k)
 230  g(j)=dum
      do 240 j=i,numscf
      dum=0.d0
      do 245 k=1,numscf
 245  dum=dum+g(k)*v(in(j)+k)
  240 h(j)=dum
      do 248 j=i,numscf
  248 v(in(i)+j)=h(j)
 250  continue
c
      cor=0.d0
      do 260 i=nap,numscf
      do 260 j=i,numscf
      x=v(in(i)+j)**2/(e(ion)+e(ion)-e(i)-e(j))
      if(j.ne.i) x=x+x
 260  cor=cor+x
      cor=cor*27.21d0
      resul=-(eiev+0.5d0*dumev+cor)
      resul2=-(eiev+dumev+cor)
      write(iw,1010)ion,eiev,dumev,cor,resul,resul2
 100  continue
      return
      end
