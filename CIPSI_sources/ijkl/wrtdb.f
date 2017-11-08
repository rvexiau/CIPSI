      subroutine wrtdb (a,n,nfile)
      implicit real*8 (a-h,o-x,z)
      dimension a(n)
      write(nfile) a
      return
      end
