      subroutine secnd(cpusec)
c     RV 01/16 rewritten with system_clock 
      integer time,ir
      real*8 cpusec
c******************************
c
c          obtain the elapsed cpu time
c
c     real*4 cpusec
cvax
c         integer*2 cptmcd,iffour,len2a(8)
c     integer*4 cputim,cptmad,izero,izero1,sys$getjpi
c     save cptmad
c     data cptmad/0/
c     equivalence (len2a(1),ifour )
c    *           ,(len2a(2),cptmcd)
c    *           ,(len2a(3),cptmad)
c    *           ,(len2a(5),izero )
c    *           ,(len2a(7),izero1)
c         data izero/0/, izero1/0/, ifour/4/, cptmcd/1031/
c      cptmad=%loc(cputim)
c      if(.not. sys$getjpi(,,,len2a,,,)) write(6,900)
c      cpusec=cptmad
c      write(6,* ) cpusec
c     cpusec=0.d0
cvax
      call system_clock(count=time, count_rate=ir)
      cpusec=real(time,kind=8)/real(ir,kind=8)
      return
900   format(22h error from sys$getjpi)
      end
