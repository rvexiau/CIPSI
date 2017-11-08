      subroutine comtab(nit,njt,iandj,tab)
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      logical*1 iandj
      dimension mmax(10),mmaz(10),tab(15,15),tac(15)
      data mmax/1,3,6,10,15,1,4,6,10,15/,mmaz/1,3,5,7,9,1,4,6,10,15/
c
      a1=2.d0/3.d0
      a2=1.d0/sqrt(5.d0)
      a3=1.5d0
      a4=sqrt(3.d0/8.d0)
      a5=sqrt(3.d0)/2.d0
      a6=0.5d0*sqrt(2.5d0)
      a7=3.d0*a2
      a8=3.d0/(2.d0*sqrt(5.d0))
c
      maxi=mmax(nit)
      maxj=mmax(njt)
      mazi=mmaz(nit)
      mazj=mmaz(njt)
c
      if(njt.le.2.or.njt.gt.5) go to 35
c     combine for j shells
      do 30 i=1,maxi
      do 10 j=1,maxj
   10 tac(j)=tab(i,j)
      if(njt.gt.3.or.njt.gt.5) go to 20
c     d shells
      tab(i,1)=tac(3)-(tac(1)+tac(2))*0.5d0
      tab(i,2)=tac(5)
      tab(i,3)=tac(6)
      tab(i,4)=(tac(1)-tac(2))*a5
      tab(i,5)=tac(4)
      go to 30
c     f shells
   20 continue
      tab(i,1)= tac(3)-(tac(5)+tac(7))*a8
      tab(i,2)=(-tac(1)-(tac(6)-4.d0*tac(8))*a2)*a4
      tab(i,3)=(-tac(2)-(tac(4)-4.d0*tac(9))*a2)*a4
      tab(i,4)=(tac(5)-tac(7))*a5
      tab(i,5)=tac(10)
      tab(i,6)=(tac(1)-tac(6)*a7)*a6
      tab(i,7)=(tac(2)-tac(4)*a7)*a6
   30 continue
   35 continue
      if(nit.le.2.or.nit.gt.5) return
      do 60 j=1,mazj
      do 40 i=1,maxi
   40 tac(i)=tab(i,j)
      if(nit.gt.3.or.nit.gt.5) go to 50
      tab(1,j)=tac(3)-(tac(1)+tac(2))*0.5d0
      tab(2,j)=tac(5)
      tab(3,j)=tac(6)
      tab(4,j)=(tac(1)-tac(2))*a5
      tab(5,j)=tac(4)
      go to 60
   50 continue
      tab(1,j)= tac(3)-(tac(5)+tac(7))*a8
      tab(2,j)=(-tac(1)-(tac(6)-4.d0*tac(8))*a2)*a4
      tab(3,j)=(-tac(2)-(tac(4)-4.d0*tac(9))*a2)*a4
      tab(4,j)=(tac(5)-tac(7))*a5
      tab(5,j)=tac(10)
      tab(6,j)=(tac(1)-tac(6)*a7)*a6
      tab(7,j)=(tac(2)-tac(4)*a7)*a6
   60 continue
      return
      end
