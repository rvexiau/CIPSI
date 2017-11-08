      function flopoi(a,n)
      implicit real*8(a-h,o-z)
      integer*2 symbl
      character*1 a(20),point,minus,plus,fig,symb
      character*1 aster,equal,blank
      logical*1 sign
      common/geod/ value(100),nsymb,symb(100,8),symbl(100)
     *,fig(16),minus,plus,equal,aster,blank
      common /opti/ dq(10),vv,g1(10),q1(10),q2(10),vv3,g3(10),vv4,g4(10)
     *,mvar(50),mvaror(50),ivar,maxopt,icompt,ico1,optiba
      point=fig(11)
      i=0
c      write(6,*)'fig'
c      write(6,21)fig,minus,plus,equal,aster,blank
 21   format(1x,80a1)
      flop=0.d0
      pow=1.d0
      sign=.false.
      if (a(1).eq.minus) i=i+1
      if (a(1).eq.minus) sign=.true.
      if (a(1).eq.plus) i=i+1
   10 i=i+1
      if (a(i).eq.point) go to 50
      do 20 j=1,10
      if (a(i).eq.fig(j)) go to 30
   20 continue
      write(6,25) (a(k),k=1,n)
   25 format (//5x,'error in function flopoi, unresolved character'/
     *        5x,50a1/)
      stop
   30 if (j.eq.10) j=0
      flop=10.d0*flop+dble(j)
      go to 10
   50 i=i+1
      if (i.gt.n) go to 80
      do 60 j=1,10
      if (a(i).eq.fig(j)) go to 70
   60 continue
      write(6,25) (a(k),k=1,n)
      stop
   70 if (j.eq.10) j=0
      pow=pow*0.1d0
      flop=flop+pow*dble(j)
      go to 50
   80 flopoi=flop
      if (sign) flopoi=-flop
      return
      end
