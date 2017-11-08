      subroutine vout(v,in,n,m)
      implicit real*8(a-h,o-z)
      dimension v(*),in(*)
      common/iofile/ir,iw
 1000 format(/)
 1001 format(5x,12(4x,i3,3x))
 1002 format(i5,12f10.6)
c
      max=12
      imax=0
   10 imin=imax+1
      imax=imax+max
      if(imax.gt.m) imax=m
      write(iw,1000)
      write(iw,1001) (i,i=imin,imax)
      write(iw,1000)
      do 20 j=1,n
   20 write(iw,1002) j,(v(in(i)+j),i=imin,imax)
      if(imax.lt.m) go to 10
      return
      end
