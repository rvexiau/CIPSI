      subroutine dnode(a,rt,k)
      implicit real*8 (a-h,o-z)
      dimension a(10),rt(10)
c
c     subroutine returns in rt(i) the ith root of a polynomial of order
c     k whose mth coefficient is stored in a(m+1).  it is assumed that
c     the initial values in rt bracket the final values.
c
      tol=1.0d-10
      k1=k+1
      r2=0.0d+00
      p2=a(1)
      do 100 m=1,k
      r1=r2
      p1=p2
      r2=rt(m)
      p2=a(k1)
      do 10 i=1,k
   10 p2=p2*r2+a(k1-i)
      prod=p1*p2
      if(prod.lt.0.0d+00) go to 20
      write(6,15) m,k
   15 format(//' root number ',i3,' was not found for polynomial of orde
     1r ',i3//)
      write(6,*) (a(i),i=1,k)
      call exit
   20 r5=r1
      p5=p1
      r6=r2
      p6=p2
   30 r3=r5
      p3=p5
      r4=r6
      p4=p6
      r =(r3*p4-r4*p3)/(p4-p3)
      dr=r4-r3
      delta=dr
      if(abs(delta).lt.tol) go to 90
      dr=0.0625d+00*dr
      r5=r-dr
      if(r5.lt.r3) r5=r3
      r6=r+dr
      if(r6.gt.r4) r6=r4
      p5=a(k1)
      p6=p5
      do 40 i=1,k
      p5=p5*r5+a(k1-i)
   40 p6=p6*r6+a(k1-i)
   45 prod=p5*p6
      if(prod.lt.0.0d+00) go to 30
      prod=p3*p5
      if(prod.gt.0.0d+00) go to 60
      r5=0.25d+00*r3+0.75d+00*r5
      p5=a(k1)
      do 50 i=1,k
   50 p5=p5*r5+a(k1-i)
      go to 45
   60 r6=0.25d+00*r4+0.75d+00*r6
      p6=a(k1)
      do 70 i=1,k
   70 p6=p6*r6+a(k1-i)
      go to 45
   90 rt(m)=r
  100 continue
      return
      end
