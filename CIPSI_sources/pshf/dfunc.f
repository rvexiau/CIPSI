      subroutine dfunc(x,n)
      implicit real*8 (a-h,o-z)
      common/ffm/ff(19)
      tol=1.0d-29
      xx=x+x
      facmin=xx
c      e=1.0d-78
      e=1.0d-38
      if(facmin.lt.360.0d+00) e=exp(-x)
      if(facmin.gt.80.0d+00) go to 100
      term=1.0d+00
      sum=1.0d+00
      fac=n
      fac=fac+0.5d+00
   10 fac=fac+1.0d+00
      term=term*x/fac
      sum=sum+term
      if(fac.le.facmin) go to 10
      t=term
      s=sum
      if(t.gt.s*tol) go to 10
      fac=n+n+1
      ff(n+1)=sum*e/fac
      m=n-1
      fac=m+m+1
   20 if(m.lt.0) return
      ff(m+1)=(e+xx*ff(m+2))/fac
      m=m-1
      fac=fac-2.0d+00
      go to 20
c
c     use asymptotic expansion for large arguments.
c
  100 a=sqrt(.7853981633974483096156608d+00/x)
      tmax=a*tol/e
      term=1.0d+00/xx
      sum=term
      fac=1.0d+00
  110 fac=fac-2.0d+00
      term=fac*term/xx
      sum=term+sum
      t=term
      if(abs(t).gt.tmax) go to 110
      ff(1)=a-e*sum
      fac=-1.0d+00
      m=0
  120 if(m.eq.n) return
      m=m+1
      fac=fac+2.0d+00
      ff(m+1)=(fac*ff(m)-e)/xx
      go to 120
      end
