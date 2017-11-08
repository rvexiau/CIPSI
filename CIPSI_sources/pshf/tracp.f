      function tracp(a,b,n)
      implicit real*8 (a-h,o-z)
      dimension a(*),b(*)
      data zero,two /0.0d+00,2.0d+00/
      tracp=zero
      m=0
      do 20 i=1,n
      do 10 j=1,i
      m=m+1
   10 tracp=tracp+a(m)*b(m)*two
   20 tracp=tracp-a(m)*b(m)
      return
      end
