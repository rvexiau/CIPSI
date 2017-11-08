      subroutine insert(c,cgel,n,nbgel,in,ind)
      implicit real*8(a-h,o-z)
      dimension ind(nbgel+1)
      dimension c(*),cgel(*),in(*)
      k=1
      jj=0
      do 10 j=1,n
          if(j.ne.ind(k)) then
              jj=jj+1
              do 20 i=1,n
                  ij=in(j)+i
                  ijj=in(jj)+i
                  c(ij)=cgel(ijj)
20            continue
           else
              k=k+1
           end if
10    continue
      return
      end
