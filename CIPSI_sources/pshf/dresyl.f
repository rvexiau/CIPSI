      subroutine dresyl (a,na,ma,b,nb,mb,n,m,test,k1)
      implicit real*8(a-h,o-z)
      common/iofile/ir,iw,ip,is
      dimension a(na,ma),b(nb,mb)
       np1=n+1
       k1=0
       do 1 l=1,n
       lp1=l+1
      if(abs(a(l,l)).gt.test) go to 2
c
c cas ou le pivot est trop petit
c
       i=l
 5     i=i+1
       if(i.gt.n) go to 101
      if(abs(a(i,l)).lt.test) go to 5
        go to 3
c
c cas ou la resolution est impossible
c
 101   k1=1
      write(iw,7)
    7 format(/' linear system for diis becomes ill-conditioned dimension
     & for extrapolation reduced ')
c
c
  11  return
   3  do 4 j=1,n
       alj=a(l,j)
      a(l,j)=a(i,j)
  4   a(i,j)=alj
       do 50 j=1,m
       blj=b(l,j)
       b(l,j)=b(i,j)
  50  b(i,j)=blj
 2     if(l.eq.n) go to 100
       do 6 j=lp1,n
   6   a(l,j)=a(l,j)/a(l,l)
 100    do 8 j=1,m
   8   b(l,j)=b(l,j)/a(l,l)
       if(l.eq.n) go to 1
       do 19 i=lp1,n
       do 70 j=lp1,n
  70  a(i,j)=a(i,j)-a(i,l)*a(l,j)
       do 80 j=1,m
  80  b(i,j)=b(i,j)-a(i,l)*b(l,j)
 19     continue
  1    continue
       i=n
 10    i=i-1
      if(i.eq.0) return
       ip1=i+1
       do 12 k=ip1,n
        do 12 j=1,m
  12  b(i,j)=b(i,j)-a(i,k)*b(k,j)
       go to 10
      end
