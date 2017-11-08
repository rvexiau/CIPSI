      subroutine focktr (h,v,f,t,ia,num)
      implicit real*8 (a-h,o-z)
      dimension h(*),v(*),f(*),t(*),ia(*)
      data zero/0.d0/
      do 70 j=1,num
      do 40 i=1,num
      dum=zero
      do 30 k=1,num
      kj=ia(j)+k
      if(k.gt.j) kj=ia(k)+j
      ik=ia(i)+k
      if(k.gt.i) ik=ia(k)+i
   30 dum=dum+f(ik)*v(kj)
   40 t(i)=dum
      do 60 i=1,j
      dum=zero
      do 50 k=1,num
      ki=ia(i)+k
      if(k.gt.i) ki=ia(k)+i
   50 dum=dum+v(ki)*t(k)
       ij=ia(j)+i
   60 h(ij)=dum
   70 continue
      return
      end
