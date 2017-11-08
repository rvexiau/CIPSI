      subroutine vecout(v,e,in,n,m)
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      character*4 iflab
      common/iofile/ir,iw
      common/runlab/iflab(3,doa)
      dimension v(*),e(*),in(*)
c
      mdim=12
      max=0
   10 min=max+1
      max=max+mdim
      if(max.gt.m) max=m
      write(iw,9999)
      write(iw,9996) (e(j),j=min,max)
      write(iw,9999)
      write(iw,9998) (j,j=min,max)
      write(iw,9999)
      do 20 i=1,n
   20 write(iw,9997) i,(iflab(k,i),k=1,3),(v(in(j)+i),j=min,max)
      if(max.lt.m) go to 10
      return
 9999 format(/)
 9998 format(16x,12(3x,i3,3x))
 9997 format(i4,3a4,12f12.5)
 9996 format(16x,12f9.4)
      end
