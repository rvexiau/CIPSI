      subroutine grossc(a,b,ia,n)
      implicit real*8 (a-h,o-z)
      dimension a(*),b(*),ia(*)
      data zero /0.0d+00/
      do 200 i=1,n
      dum=zero
      do 100 j=1,n
      ij=ia(i)+j
      if(j.gt.i) ij=ia(j)+i
  100 dum=dum+a(ij)
  200 b(i)=dum
      return
      end
