      subroutine droot
c     this version uses christoffel formula for weights.
      implicit real*8 (a-h,o-z)
      real*8 xx,uf,wf
      common/root/xx,uf(9),wf(9),nroots
      common/ffm/ff(19)
      common/roots/r(9,9),w(9,9)
      dimension c(10,10),s(10,10),a(10),rt(10)
c
c     ith root of the jth rys polynomial is returned in r(i,j) with
c     the corresponding weight factor in w(i,j).   j=1,2,...,n
c
      n=nroots
      x=xx
      if(n.lt.2) n=2
      n1=n+1
      nn=n+n
      call dfunc(x,nn)
      do 10 i=1,n1
      do 10 j=1,n1
   10 s(i,j)=ff(i+j-1)
      call dsmit(c,s,n1)
      do 20 i=1,n
      do 20 j=1,i
      w(i,j)=0.0d+00
   20 r(i,j)=0.0d+00
      wsum=ff(1)
      w(1,1)=wsum
      r(1,1)=ff(2)/wsum
      dum=sqrt(c(2,3)**2-4.0d+00*c(1,3)*c(3,3))
      r(1,2)=0.5d+00*(-c(2,3)-dum)/c(3,3)
      r(2,2)=0.5d+00*(-c(2,3)+dum)/c(3,3)
      if(n.eq.2) go to 70
      do 25 i=3,n1
   25 rt(i)=1.0d+00
      rt(1)=r(1,2)
      rt(2)=r(2,2)
      do 60 k=3,n
      k1=k+1
      do 30 i=1,k1
   30 a(i)=c(i,k1)
c      write(6,*) 'droot', (a(i),i=1,k1)
      call dnode(a,rt,k)
      do 50 i=1,k
   50 r(i,k)=rt(i)
   60 continue
   70 do 150 k=2,n
      jmax=k-1
      do 150 i=1,k
      root=r(i,k)
      dum=1.0d+00/ff(1)
      do 110 j=1,jmax
      j1=j+1
      poly=c(j1,j1)
      do 100 m=1,j
  100 poly=poly*root+c(j1-m,j1)
  110 dum=dum+poly*poly
  150 w(i,k)=1.d+00/dum
      do 160 k=1,nroots
      dum=r(k,nroots)
      uf(k)=dum/(1.0d+00-dum)
  160 wf(k)=w(k,nroots)
      return
      end
