      subroutine elimat(f,nbgel,ind,n,ia,iga)
      implicit real*8(a-h,o-z)
      include 'pshf.prm'
      dimension ind(nbgel+1),ia(*),iga(*)
      dimension f(*)
      dimension temp(doas)
      jj=0
      kj=1
      do 10 j=1,n
          if(j.ne.ind(kj)) then
              ii=0
              jj=jj+1
              ki=1
              do 20 i =1,j
                  if(i.ne.ind(ki)) then
                      ii=ii+1
                      iijj=iga(jj)+ii
                      if(ii.gt.jj) iijj=iga(ii)+jj
                      ij=ia(j)+i
                      temp(iijj)=f(ij)
                   else
                      ki=ki+1
                   end if
20            continue
           else
              kj=kj+1
           end if
10    continue
      do 30 j=1,n-nbgel
          do 40 i=1,j
              ij=iga(j)+i
              f(ij)=temp(ij)
40        continue
30    continue
      return
      end
