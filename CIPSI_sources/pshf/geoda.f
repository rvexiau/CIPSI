      subroutine geoda
      implicit real*8 (a-h,o-z)
      character*1 a(80),aa(20),symb,name(100,8)
      character*1 aster,blank,minus,plus,equal,figure
      integer*2 symbl,namel(100)
      logical*1 deter,sign,rx,ry,rz
      dimension title(10)
      common/geod/ value(100),nsymb,symb(100,8),symbl(100)
     *,figure(16),minus,plus,equal,aster,blank
      common /opti/ dq(10),vv,g1(10),q1(10),q2(10),vv3,g3(10),vv4,g4(10)
     *,ivar,maxopt,icompt,ico1,optiba
      common/opti7/ww,fin
      logical*1 ww,fin
      common/iofile/ ir,iw,ip,is,iq,ih,iv
      common /zmat/ iz(50,4),bl(50),alpha(50),beta(50),nz,nzz,
     * ipar(15,5),nipar(5),npar,nstep,xx(5),dx(5),num,numb,confac,
     * name,namel,rx,ry,rz
      data a/'1','2','3','4','5','6','7','8','9','0','.',69*'0'/
      data aa /'-','+','=','*',' ',15*'0'/
      write(6,*)'entree dans geoda'
      do 1 i=1,11
 1    figure(i)=a(i)
      minus=aa(1)
      plus=aa(2)
      equal=aa(3)
      aster=aa(4)
      blank=aa(5)
      confac=1.d0
      rx=.false.
      ry=.false.
      rz=.false.
      ww=.true.
c      open (unit=53,status='scratch')
      ib=53
c
c     skip z-matrix
c
   10 read (5,20) a
c      write(6,21) a
      write(ib) a
   20 format (80a1)
   21 format (1x,80a1)
      ia=0
   30 ia=ia+1
      if (ia.gt.80) go to 50
      if (a(ia).eq.blank) go to 30
      go to 10
c
c     read the values of the parameters
c
   50 nsymb=0
   60 read (5,20,end=100) a
c      write(6,*)'val des param'
c      write(6,21) a
c
c     ... skip blanks
c
      ia=0
   70 ia=ia+1
      if (ia.gt.80) go to 130
      if (a(ia).eq.blank) go to 70
c
c     ... store symbolic name
c
   80 nsymb=nsymb+1
      j=0
   90 j=j+1
      symb(nsymb,j)=a(ia)
c      write(6,*)' symb name'
c      write(6,21) a(ia)
      ia=ia+1
      if (a(ia).ne.blank.and.a(ia).ne.equal) go to 90
      symbl(nsymb)=j
c
c     ... skip blanks and = sign
c
  100 ia=ia+1
      if (a(ia).eq.blank.or.a(ia).eq.equal) go to 100
c
c     ... store value
c
      j=0
  120 j=j+1
      aa(j)=a(ia)
      ia=ia+1
c      write(6,*) 'aa(j)'
c      write(6,21) aa(j)
      if (a(ia).ne.blank) go to 120
c      write(6,21) aa
      value(nsymb)=flopoi(aa,j)
      go to 60
  130 continue
      do 132 i=1,nsymb
      lng=symbl(i)
      nbl=9-lng
  132 if(ww)write(6,135) (symb(i,j),j=1,lng),(blank,j=1,nbl),value(i)
  135 format (/5x,9a1,'=',f20.12)
c
c     read z matrix
c
      rewind ib
      iz(1,1)=0
      iz(1,2)=0
      iz(1,3)=0
      iz(1,4)=0
      bl(1)=0.d0
      alpha(1)=0.d0
      beta(1)=0.d0
      nz=0
 150  read(ib) a
c      write(6,21) a
c
c     ... skip blanks
c
      ia=0
  160 ia=ia+1
      if (ia.gt.80) go to 500
      if (a(ia).eq.blank) go to 160
c
c     ... store centre name
c
      nz=nz+1
      j=0
  170 j=j+1
      name(nz,j)=a(ia)
      ia=ia+1
      if (a(ia).ne.blank) go to 170
      namel(nz)=j
      if (nz.eq.1) go to 150
c
c     ... skip blanks
c
  180 ia=ia+1
      if (a(ia).eq.blank) go to 180
c
c     ... find iz(1)
c
      j=0
  190 j=j+1
      aa(j)=a(ia)
      ia=ia+1
      if (a(ia).ne.blank) go to 190
      call find(aa,j,name,namel,nz-1,n)
      iz(nz,1)=n
c
c     ... skip blanks
c
  200 ia=ia+1
      if (a(ia).eq.blank) go to 200
c
c     ... find the value of bl
c
      sign=.false.
      if (a(ia).ne.minus) go to 210
      sign=.true.
      ia=ia+1
  210 if (a(ia).eq.plus) ia=ia+1
      deter=.false.
      do 220 i=1,11
  220 if (a(ia).eq.figure(i)) deter=.true.
      j=0
  230 j=j+1
      aa(j)=a(ia)
      ia=ia+1
      if (a(ia).ne.blank) go to 230
      if (deter) bl(nz)=flopoi(aa,j)
      if (.not.deter) call find(aa,j,symb,symbl,nsymb,n)
      if (.not.deter) bl(nz)=value(n)
      if (sign) bl(nz)=-bl(nz)
      if (nz.eq.2) go to 150
c
c     ... skip blanks
c
  240 ia=ia+1
      if (a(ia).eq.blank) go to 240
c
c     ... find iz(2)
c
      j=0
  250 j=j+1
      aa(j)=a(ia)
      ia=ia+1
      if (a(ia).ne.blank) go to 250
      call find(aa,j,name,namel,nz-1,n)
      iz(nz,2)=n
c
c     ... skip blanks
c
  260 ia=ia+1
      if (a(ia).eq.blank) go to 260
c
c     ... find the value of alpha
c
      sign=.false.
      if (a(ia).ne.minus) go to 270
      sign=.true.
      ia=ia+1
  270 if (a(ia).eq.plus) ia=ia+1
      deter=.false.
      do 280 i=1,11
  280 if (a(ia).eq.figure(i)) deter=.true.
      j=0
  290 j=j+1
      aa(j)=a(ia)
      ia=ia+1
      if (a(ia).ne.blank) go to 290
      if (deter) alpha(nz)=flopoi(aa,j)
      if (.not.deter) call find(aa,j,symb,symbl,nsymb,n)
      if (.not.deter) alpha(nz)=value(n)
      if (sign) alpha(nz)=-alpha(nz)
      if (nz.eq.3) go to 150
c
c     ... skip blanks
c
  300 ia=ia+1
      if (a(ia).eq.blank) go to 300
c
c     ... find iz(3)
c
      j=0
  310 j=j+1
      aa(j)=a(ia)
      ia=ia+1
      if (a(ia).ne.blank) go to 310
      call find(aa,j,name,namel,nz-1,n)
      iz(nz,3)=n
c
c     ... skip blanks
c
  320 ia=ia+1
      if (a(ia).eq.blank) go to 320
c
c     ... find the value of beta
c
      sign=.false.
      if (a(ia).ne.minus) go to 330
      sign=.true.
      ia=ia+1
  330 if (a(ia).eq.plus) ia=ia+1
      deter=.false.
      do 340 i=1,11
  340 if (a(ia).eq.figure(i)) deter=.true.
      j=0
  350 j=j+1
      aa(j)=a(ia)
      ia=ia+1
      if (a(ia).ne.blank) go to 350
      if (deter) beta(nz)=flopoi(aa,j)
      if (.not.deter) call find(aa,j,symb,symbl,nsymb,n)
      if (.not.deter) beta(nz)=value(n)
      if (sign) beta(nz)=-beta(nz)
c
c     ... skip last blanks and set iz(4)
c
      iz(nz,4)=0
  360 ia=ia+1
      if (ia.gt.80) go to 150
      if (a(ia).eq.blank) go to 360
      iz(nz,4)=1
      if (a(ia).eq.minus) iz(nz,4)=-1
      go to 150
c
c     end of z-matrix
c
  500 continue
      call buildz
      return
      end
