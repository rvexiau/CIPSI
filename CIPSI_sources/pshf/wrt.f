      subroutine wrt(index)
      implicit real*8 (a-h,o-z)
      common/iofile/ir,iw
      common/wrtc/i1,i2,i3,i4,inx(4,5),q4,val,q4x(5),valx(5),nwrt
 9999 format(1x,5(4i3,f5.3,f9.6))
      if(index.eq.1) go to 10
      nwrt=nwrt+1
      inx(1,nwrt)=i1
      inx(2,nwrt)=i2
      inx(3,nwrt)=i3
      inx(4,nwrt)=i4
      q4x(nwrt)=q4
      valx(nwrt)=val
      if(nwrt.lt.5) return
   10 write(iw,9999) ((inx(i,j),i=1,4),q4x(j),valx(j),j=1,nwrt)
      nwrt=0
      return
      end
