      subroutine getgel(cc,v,n,nbgel,ind,in)
      implicit real*8(a-h,o-z)
      dimension cc(*),v(*),ind(nbgel+1),in(*)
      do 10 j=1,nbgel
          do 20 i=1,n
             ij=in(j)+i
             ijj=in(ind(j))+i
             v(ij)=cc(ijj)
20         continue
10    continue
      return
      end
