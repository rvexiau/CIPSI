      subroutine readb (a,n,nfile)
      implicit real*8 (a-h,o-x,z)
      dimension a(n)
      read(nfile) a
      return
      end
