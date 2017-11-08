      subroutine dmtx(v,f,q,num)
      implicit real*8 (a-h,o-z)
      dimension v(*),f(*),q(*)
      data zero/0.d0/
      ni=-num
      jk=0
      do 5 j=1,num
      do 5 k=1,j
      jk=jk+1
    5 f(jk)=zero
      do 10 i=1,num
      ni=ni+num
      if(q(i).eq.zero) go to 10
      jk=0
      do 20 j=1,num
      do 20 k=1,j
      jk=jk+1
   20 f(jk)=f(jk)+v(ni+j)*v(ni+k)*q(i)
   10 continue
      return
      end
