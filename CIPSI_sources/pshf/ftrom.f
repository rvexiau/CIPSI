      subroutine ftrom(c,f,tf,n,in,ia)
      implicit real*8(a-h,o-z)
      include 'pshf.prm'
      dimension c(*),f(*),in(*),ia(*),tf(*),t(doa)
      do 10 j=1,n
          do 20 k=1,n
             dum=0.d0
             do 30 l=1,n
                kl=ia(l)+k
                if(k.gt.l) kl=ia(k)+l
                lj=in(j)+l
                dum=dum+f(kl)*c(lj)
30           continue
             t(k)=dum
20         continue
           do 40 i=1,n
               dum=0.d0
               do 50 k=1,n
                   ki=in(i)+k
                   dum=dum+c(ki)*t(k)
50             continue
               ij=ia(j)+i
               if(i.gt.j) ij=ia(i)+j
               tf(ij)=dum
40         continue
10    continue
      return
      end
