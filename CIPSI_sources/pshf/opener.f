      subroutine opener
      implicit real*8 (a-h,o-z)
      dimension enx(300)
      rewind 33
      do i=1,300
      enx(i)=0.
      end do
      write(33) enx
      rewind 33
      write (6,*) ' file 33 initialisee a zero'
      return
      end
