      subroutine print(n1,n2)
      implicit real*8 (a-h,o-z)
      common/iofile/ir,iw,ip,is
      dimension nn(48)
      common/symmat/t(432)
      imax=n1-1
  100 imin=imax+1
      imax=imax+4
      if(imax.gt.n2) imax=n2
      nj=9*n1-8
      do 200 j=1,3
      ni=0
      do 150 i=imin,imax
      nn(i)=nj+ni
  150 ni=ni+9
      write(iw,1000) (t(nn(i)),t(nn(i)+1),t(nn(i)+2),
     1 i=imin,imax)
  200 nj=nj+3
      write(iw,1001)
      if(imax.lt.n2) go to 100
 1000 format(4x,4(3f10.5,2h *))
 1001 format(/)
      return
      end
