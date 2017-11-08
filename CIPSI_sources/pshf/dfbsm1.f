      subroutine dfbsm1
      implicit real*8(a-h,o-z)
      common/uncp/ax,ay,az,bx,by,bz,xc,yc,zc,ga,gb,gc,am,bm,cm,
     1 pt,alfa,alfb,alfc,ctta,w,f1,f2,f3,f4,f5,f6,f7,f8,f9,xbs,bs(10),
     2 xg(64),hg(64),rxg(40),h1(40),h23(40),h456(40),h78(40),h9(40),
     3 fgam(13),rab,idum(4),nbs
      dimension dntab(10),an2tab(10)
      data dntab/0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0/
      data an2tab/2.0,6.0,10.,14.,18.,22.,26.,30.,34.,38./
c...... fonction de bessel spherique modifiee de premiere espece
c...... parametres nbs,xbs
c...... resultat dans bs(nbs+1)
      n1=nbs+1
      dn=dntab(n1)
      an2=an2tab(n1)
      dabx=xbs
      if(xbs.lt.0.0) dabx=-xbs
      if(dabx.gt.an2) go to 30
c......
      t=9.45d0*dabx
      an=dsqrt(t)-0.5d0
      p=dn
      if(dn.lt.an) p=an
      m=t/(p+0.5d0)+p+3.7d0
      r=0.0d0
      do 10 j=1,m
      k=m-j+2
      ak=k+k-1
      r=xbs/(xbs*r+ak)
      if(n1.lt.k) go to 10
      bs(k)=r
   10 continue
      bs(1)=1.0d0/(r*xbs+dabx+1.0d0)
      if(nbs.eq.0) return
      do 20 j=2,n1
      bs(j)=bs(j-1)*bs(j)
   20 continue
      return
c......
   30 tx=dabx+dabx
      te=0.0d0
      if(tx.lt.87.3d0) te=dexp(-tx)
      bs(1)=(1.0d0-te)/tx
      if(nbs.eq.0) return
      bs(2)=(0.5d0*te+0.5d0-bs(1))/xbs
      if(nbs.eq.1) return
      do 40 j=3,n1
      ajm3=j+j-3
   40 bs(j)=bs(j-2)-ajm3*bs(j-1)/xbs
      return
      end
