      subroutine vec(u,c,j,k)
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      dimension c(3,dc),r(3),u(3)
      zero=0.0d0
      r2=zero
      do 10 i=1,3
      r(i)=c(i,j)-c(i,k)
   10 r2=r2+r(i)*r(i)
      r2=sqrt(r2)
      do 20 i=1,3
   20 u(i)=r(i)/r2
      return
      end
