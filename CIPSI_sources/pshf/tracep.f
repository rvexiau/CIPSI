      function tracep(a,b,n)
      implicit real*8(a-h,o-z)
      dimension a(*),b(*)
      data zero,two/0.d0,2.d0/
      tracep=zero
      k=0
      do 20 i=1,n
      do 10 j=1,i
      k=k+1
   10 tracep=tracep+a(k)*b(k)*two
   20 tracep=tracep-a(k)*b(k)
      return
      end
