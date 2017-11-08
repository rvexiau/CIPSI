      subroutine dsmit(c,s,n)
      implicit real*8 (a-h,o-z)
      dimension c(10,10),s(10,10),v(10),y(10)
c
c     subroutine returns an n by n triangular matrix c such that
c     c(transpose)sc=i,  where i is an n by n identity matrix.
c
      do 10 i=1,n
      do 10 j=1,i
   10 c(i,j)=0.0d+00
      do 100 j=1,n
      kmax=j-1
      fac=s(j,j)
      if(kmax.eq.0) go to 60
      do 20 k=1,kmax
      v(k)=0.0d+00
   20 y(k)=s(k,j)
      do 50 k=1,kmax
      dot=0.0d+00
      do 30 i=1,k
   30 dot=c(i,k)*y(i)+dot
      do 40 i=1,k
   40 v(i)=v(i)-dot*c(i,k)
   50 fac=fac-dot*dot
   60 fac=1.0d+00/sqrt(fac)
      c(j,j)=fac
      if(kmax.eq.0) go to 100
      do 70 k=1,kmax
   70 c(k,j)=fac*v(k)
  100 continue
      return
      end
