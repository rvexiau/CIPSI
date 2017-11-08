      subroutine elimv(v,velim,nbgel,ind,in,n)
      implicit real*8(a-h,o-z)
      dimension v(*),velim(*),in(*)
      dimension ind(nbgel+1)
      k=1
      jj=0
      do 10 j=1,n
         if(j.ne.ind(k)) then
            jj=jj+1
            do 20 i=1,n
                ijj=in(jj)+i
                ij=in(j)+i
                velim(ijj)=v(ij)
20          continue
         else
            k=k+1
         end if
10    continue
      return
      end
