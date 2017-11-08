      subroutine addmat(a,b,c,n)
      implicit real*8(a-h,o-z)
      dimension a(*),b(*),c(*)
      do 10 i=1,n
   10 c(i)=a(i)+b(i)
      return
      end
