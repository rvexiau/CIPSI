      program rcut
Ccc   RG : modified on 25/11/08 : if the cutoff radii are negative, now calculates
Ccc   RG : modified on 25/11/08 : a difference with respect to a rpol calculation
      implicit real*8(a-h,o-z)
      integer plder,plmax
      parameter (plder=3)
      include 'pshf.prm'
      logical diffCalc,ysoft,firstRun
      common/rayona/rcut1(dc,0:plder)
      common/puisno/no(2)
      common /diff/ diffCalc,firstRun
      namelist/rcutval/plmax,rcut1,ysoft,no
      character *26 timeday
      call openf
      call fdate(timeday)
      write(6,9999)timeday
      diffCalc=.false.
      firstRun=.false.
      read(*,rcutval)
      rewind 5
      if(rcut1(1,0)<0.d0) then
        diffCalc=.true.
      endif
cCc ============================
cCc = Differential calculation =
cCc ============================      
      if(diffCalc) then
        write(*,*) '********** Differential calculation **********'
        write(*,*) '**********************************************'
        firstRun=.true.
        call set
        call gen
        call efield
        firstRun=.false.
      endif
      call set
      call gen
      call efield
      call fdate(timeday)
      write(6,9998)timeday
 9999 format(/,80(1h*),/,8(10h   rcut   ),/, 10x,a10,5x,a10,/,
     * 80(1h*))
 9998 format(/,80(1h*),/,8(10h fin rcut ),/, 10x,a10,5x,a10,/,
     * 80(1h*))
      stop
      end
