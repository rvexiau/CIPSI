      subroutine polzer
      implicit real*8(a-h,o-z)
      include 'pshf.prm'
      character*40 typener
      logical ymono
      common/iofile/ir,iw
      common/polar/calfa(dc)
      common/polxx/polz
      common/infoa/nat,ich,mul,num,nnp,ne,na,nb,zan(dc),c(3,dc),ymono
      dimension el(3)
      common/elchp/elec(3)
c
c    elec   champ constant additionel pour
c
      polz=0.0d0
      do 1000 ic=1,nat
      alfad=calfa(ic)*0.5d0
      do 5 k=1,3
    5 el(k)=elec(k)
      do 50 i=1,nat
      if(i.eq.ic) go to 50
      ri=0.d0
      do 10 k=1,3
   10 ri=ri+(c(k,ic)-c(k,i))**2
      ri=ri**1.5d0
      ri=1.d0/ri
      zi=zan(i)
      do 20 k=1,3
c champ cree sur l'atome ic par les autres atomes
   20 el(k)=el(k)+zi*(c(k,ic)-c(k,i))*ri
   50 continue
      dum=0.d0
      do 60 k=1,3
   60 dum=dum+el(k)**2
      epz=-dum*alfad
      write(iw,9999) ic,(el(k),k=1,3),epz
      polz=polz+epz
 1000 continue
      write(iw,9998) polz
 9999 format(/5x,38h nuclei core for polarization for atom,i3,5x,
     1 4d15.8)
 9998 format(/5x,31h total nuclei core polarization,/,10x,d15.8)
      typener='energie de polarisation des coeurs'
      call stkener (typener,polz,1)
      return
      end
