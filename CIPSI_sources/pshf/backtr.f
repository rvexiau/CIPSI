      subroutine backtr(v,f,t,ia,in,num)
      implicit real*8(a-h,o-z)
      dimension v(*),f(*),t(*),ia(*),in(*)
      data zero/0.d0/
      do 30 j=1,num
      do 20 i=1,num
      dum=zero
c
      do 10 k=1,num
      jk=in(j)+k
      if(k.gt.i) go to 5
      ki=ia(i)+k
      go to 10
    5 ki=ia(k)+i
   10 dum=dum+f(ki)*v(jk)
   20 t(i)=dum
      do 30 i=1,num
   30 v(in(j)+i)=t(i)
      return
      end
