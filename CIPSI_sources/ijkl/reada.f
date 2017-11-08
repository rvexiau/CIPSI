      subroutine reada (a,n,nav)
      implicit real*8 (a-h,o-x,z)
      dimension a(n)
      read(10,rec=nav) a
      return
      end
