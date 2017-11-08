      subroutine hdiag(c,t,e,ia,in,no)
      implicit real*8(a-h,o-z)
      dimension c(*),t(*),ia(*),in(*),e(*)
c     methode de jacobi
c      calcul des valeurs propres et des vecteurs propres
      write(*,*) 'dimensions in hdiag = ',no,in(no)+no,ia(no)+no
      n=no
      nx=n*(n+1)/2
      r2=sqrt(.5d0)
      delta=0.1d0
      do 10 i=1,n
      do 10 j=1,n
      t(in(i)+j)=0.d0
      if(i.eq.j) t(in(i)+j)=1.d0
 10   continue
 35   n1=n-1
      i4=0
      do 85 i=1,n1
      i2=i+1
      do 85 j=i2,n
      a7=c(ia(j)+i)
      a9=a7*a7
      if(a9-delta)85,85,40
 40   i4=1
      a1=c(ia(i)+i)
      a2=c(ia(j)+j)
      a3=a2-a1
      if(a3)50,45,50
 45   a3=r2
      a4=a3
      a9=.5d0
      a10=a9
      a6=a7
      goto57
 50   a7=2.d0*a7
      et=a7/a3
      a=et/(sqrt(1.d0+et*et)+1.d0)
      a9=1.d0/(1.d0+a*a)
      a10=1.d0-a9
      a3=sqrt(a9)
      a4=a*a3
      a6=a7*a*a9
 57   a5=a1*a9+a2*a10
      a8=a1*a10+a2*a9
      do 84 k=1,n
      if(i-k)60,70,60
 60   if(j-k)65,75,65
 65   continue
      ik=ia(i)+k
      if(k.gt.i) ik=ia(k)+i
      jk=ia(j)+k
      if(k.gt.j) jk=ia(k)+j
      a9=c(ik)
      a14=c(jk)
      c(ik)=a9*a3-a14*a4
      c(jk)=a9*a4+a14*a3
      goto 80
 70   c(ia(i)+k)=a5-a6
      c(ia(j)+i)=0.d0
      goto 80
 75   c(ia(j)+k)=a8+a6
 80   a15=t(in(i)+k)
      a16=t(in(j)+k)
      t(in(i)+k)=a15*a3-a16*a4
      t(in(j)+k)=a15*a4+a16*a3
 84   continue
 85   continue
      if(i4)35,90,35
 90   if(delta-1.d-18)100,100,95
 95   continue
      delta=delta/10.
      goto 35
 100  continue
      do 101 i=1,n
  101 e(i)=c(ia(i)+i)
  102 continue
      n1=n-1
      do 115 i=1,n1
      if(e(i)-e(i+1))115,115,105
 105  a=e(i)
      e(i)=e(i+1)
      e(i+1)=a
      do 110 k=1,n
      a=t(in(i)+k)
      t(in(i)+k)=t(in(i+1)+k)
 110  t(in(i+1)+k)=a
      goto 102
 115  continue
      return
      end
