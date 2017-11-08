      subroutine timit(index)
      implicit real*8 (a-h,o-z)
      common/iofile/ir,iw,ip,is
      common/times/ti,tx,tim,tom,to
      call secnd(tom)
      tim=tom-to
      tx=tim-ti
      ti=tim
      if(index.ne.0) write(iw,1000) tx,tim
 1000 format(10x,'elapsed time = ',f10.3,5x,'total time = ',f10.3)
      return
      end
