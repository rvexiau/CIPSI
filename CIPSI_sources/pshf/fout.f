      subroutine fout(a,n)
      implicit real*8(a-h,o-z)
      common/iofile/ir,iw
      dimension a(*)
 1000 format(/)
 1001 format(4x,14(3x,i3,3x))
 1002 format(i4,14f9.5)
      max=14
      imax=0
  100 imin=imax+1
      imax=imax+max
      if(imax.gt.n) imax=n
      write(iw,1000)
      write(iw,1001) (i,i=imin,imax)
      write(iw,1000)
      do 400 j=imin,n
      isup=imax
      if(isup.gt.j) isup=j
      ni=j*(j-1)/2
      write(iw,1002) j,(a(ni+i),i=imin,isup)
  400 continue
      if(imax.lt.n) go to 100
      return
      end
