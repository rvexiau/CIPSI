      subroutine combin(ns1,n,m,nt,md,nd)
      integer*4 ns1(md,nd)
      i=1
      do 1 j=1,m
1      ns1(i,j)=j
5     continue
      j1=0
      do 3 j=1,m
 3     if(ns1(i,j).eq.n-m+j) j1=j1+1
      if(j1.eq.m) goto 4
      i=i+1
      do 2 j=1,m
2     ns1(i,j)=ns1(i-1,j)
      ns1(i,m-j1)=ns1(i-1,m-j1)+1
      if(j1.eq.0) goto 5
      mmj1=m-j1+1
      do 6 j=mmj1,m
6     ns1(i,j)=ns1(i,j-1)+1
      goto 5
4      continue
      nt=i
      return
      end
