      subroutine reduc
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
c     this subroutine reduces d and f transformation tables from 6 0r 10
c     components to 5 or 7
      common/symtry/ invt(48),iso(ds,12),ptr(3,144),dtr(6,288),
     1 ftr(10,480),nt
      dimension d5d6(6,6),d6d5(6,6),f7f10(10,10),f10f7(10,10)
      dimension dx(6,6),fx(10,10)
      equivalence (dx(1,1),fx(1,1))
c
      dsq2=sqrt(2.d0)
      dsq3=sqrt(3.d0)
      dsq5=sqrt(5.d0)
      dsq7=sqrt(7.d0)
      dsq10=sqrt(10.d0)
c
      do 1 i=1,6
      do 1 j=1,6
      d5d6(i,j)=0.d0
    1 d6d5(i,j)=0.d0
      do 2 i=1,10
      do 2 j=1,10
      f7f10(i,j)=0.d0
    2 f10f7(i,j)=0.d0
c
      d5d6(1,4)=dsq3*0.5d0
      d5d6(2,4)=-d5d6(1,4)
      d5d6(1,1)=-0.5d0
      d5d6(2,1)=-0.5d0
      d5d6(3,1)=1.d0
      d5d6(5,2)=1.d0
      d5d6(6,3)=1.d0
      d5d6(4,5)=1.d0
      d5d6(1,6)=1.d0/dsq5
      d5d6(2,6)=d5d6(1,6)
      d5d6(3,6)=d5d6(1,6)
c
      d6d5(4,1)=1.d0/dsq3
      d6d5(1,1)=-1.d0/3.d0
      d6d5(6,1)=dsq5/3.d0
      d6d5(4,2)=-d6d5(4,1)
      d6d5(1,2)=d6d5(1,1)
      d6d5(6,2)=d6d5(6,1)
      d6d5(1,3)=2.d0/3.d0
      d6d5(6,3)=d6d5(6,1)
      d6d5(5,4)=1.d0
      d6d5(2,5)=1.d0
      d6d5(3,6)=1.d0
c
      f7f10(3,1)=1.d0
      f7f10(5,1)=-1.5d0/dsq5
      f7f10(7,1)=-1.5d0/dsq5
      f7f10(1,2)=-dsq3*0.5d0/dsq2
      f7f10(6,2)=-dsq3*0.5d0/dsq10
      f7f10(8,2)=dsq3*2.d0/dsq10
      f7f10(2,3)=f7f10(1,2)
      f7f10(4,3)=f7f10(6,2)
      f7f10(9,3)=f7f10(8,2)
      f7f10(5,4)=dsq3*0.5d0
      f7f10(7,4)=-f7f10(5,4)
      f7f10(10,5)=1.d0
      f7f10(1,6)=dsq5*0.5d0/dsq2
      f7f10(6,6)=-1.5d0/dsq2
      f7f10(2,7)=f7f10(1,6)
      f7f10(4,7)=f7f10(6,6)
      f7f10(1,8)=dsq5/dsq7
      f7f10(6,8)=1.d0/dsq7
      f7f10(8,8)=f7f10(6,8)
      f7f10(2,9)=f7f10(1,8)
      f7f10(4,9)=f7f10(6,8)
      f7f10(9,9)=f7f10(6,8)
      f7f10(3,10)=f7f10(1,8)
      f7f10(5,10)=f7f10(6,8)
      f7f10(7,10)=f7f10(6,8)
c
      f10f7(2,1)=-dsq2*dsq3*0.1d0
      f10f7(6,1)=1.d0/dsq10
      f10f7(8,1)=0.6d0*dsq7/dsq5
      f10f7(3,2)=f10f7(2,1)
      f10f7(7,2)=f10f7(6,1)
      f10f7(9,2)=f10f7(8,1)
      f10f7(1,3)=2.d0/5.d0
      f10f7(10,3)=f10f7(8,1)
      f10f7(3,4)=-dsq2*0.5d0/(dsq3*dsq5)
      f10f7(7,4)=-1.d0/dsq2
      f10f7(9,4)=dsq7*0.2d0
      f10f7(1,5)=-1.d0/dsq5
      f10f7(4,5)=1.d0/dsq3
      f10f7(10,5)=dsq7*0.2d0
      f10f7(2,6)=f10f7(3,4)
      f10f7(6,6)=f10f7(7,4)
      f10f7(8,6)=f10f7(9,4)
      f10f7(1,7)=f10f7(1,5)
      f10f7(4,7)=-f10f7(4,5)
      f10f7(10,7)=f10f7(10,5)
      f10f7(2,8)=2.d0*dsq2/(dsq3*dsq5)
      f10f7(8,8)=dsq7*0.2d0
      f10f7(3,9)=f10f7(2,8)
      f10f7(9,9)=f10f7(8,8)
      f10f7(5,10)=1.d0
c
      do 100 it=1,nt
      nd=6*(it-1)
      nf=10*(it-1)
c     write(6,*) it
c1020 format(10f8.3)
      do 10 i=1,6
      do 10 j=1,6
      dx(i,j)=0.d0
      do 10 k=1,6
   10 dx(i,j)=dx(i,j)+dtr(i,nd+k)*d5d6(k,j)
      do 20 i=1,6
      do 20 j=1,6
      dtr(i,nd+j)=0.d0
      do 20 k=1,6
   20 dtr(i,nd+j)=dtr(i,nd+j)+d6d5(i,k)*dx(k,j)
c
      do 30 i=1,10
      do 30 j=1,10
      fx(i,j)=0.d0
      do 30 k=1,10
   30 fx(i,j)=fx(i,j)+ftr(i,nf+k)*f7f10(k,j)
      do 40 i=1,10
      do 40 j=1,10
      ftr(i,nf+j)=0.d0
      do 40 k=1,10
   40 ftr(i,nf+j)=ftr(i,nf+j)+f10f7(i,k)*fx(k,j)
c     do 50 i=1,10
c  50 write(6,1020) (ftr(i,nf+j),j=1,10)
  100 continue
      return
      end
