      subroutine spin (sz,s2)
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
c
c     ----- calculate expectation value of -sz- and -s**2
c           for the unrestricted hf wavefunction -----
c
      common/iofile/ir,iw,ip,is,iq,ih,iv,if2,if3,if4,ioda(19)
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(dc),c(3,dc)
      common/hfden/da(doas),db(doas),ia(doa)
      dimension s(doas),t(doa)
      data zero,two /0.0d+00,2.0d+00/
 9999 format(/,10x,19(1h-),/,10x,'spin sz   = ',f7.5,/,
     1 10x,'s-squared = ',f7.5,/,10x,19(1h-))
c
      do 10 i=1,num
   10 ia(i)=i*(i-1)/2
c     ----- read in overlap and density matrices (alpha+beta) -----
c
      call reada(s,nx,2)
      call reada(da,nx,6)
c
c     ----- d = s*da*s -----
c
      do 1000 j=1,num
      do 600 i=1,num
      dum=zero
      do 500 k=1,num
      if(k.gt.i) go to 100
      ik=ia(i)+k
      go to 200
  100 ik=ia(k)+i
  200 if(k.gt.j) go to 300
      jk=ia(j)+k
      go to 400
  300 jk=ia(k)+j
  400 dum=dum+da(ik)*s(jk)
  500 continue
  600 t(i)=dum
      do 1000 i=1,j
      dum=zero
      do 900 k=1,num
      if(k.gt.i) go to 700
      ik=ia(i)+k
      go to 800
  700 ik=ia(k)+i
  800 dum=dum+s(ik)*t(k)
  900 continue
      ij=ia(j)+i
 1000 db(ij)=dum
c
c
      call reada(da,nx,9)
c
c     ----- calculate spin quantum numbers -----
c
      sz=dble(na-nb)/two
      s2=sz*sz+dble(na+nb)/two-tracep(da,db,num)
      write(iw,9999) sz,s2
      return
      end
