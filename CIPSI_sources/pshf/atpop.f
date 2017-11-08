      subroutine atpop(a,ia,b,nat)
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      common/atlim/limlow(dc),limsup(dc)
      dimension a(*),b(*),ia(*)
      data zero /0.0d+00/
      do 200 i=1,nat
      do 200 j=1,nat
      i1=limlow(i)
      i2=limsup(i)
      j1=limlow(j)
      j2=limsup(j)
      dum=zero
      if(i1.eq.0.or.j1.eq.0)go to 110
      do 100 k=i1,i2
      do 100 l=j1,j2
      kl=ia(k)+l
      if(l.gt.k) kl=ia(l)+k
  100 dum=dum+a(kl)
  110 continue
      ij=ia(i)+j
      if(j.gt.i) ij=ia(j)+i
  200 b(ij)=dum
      return
      end
