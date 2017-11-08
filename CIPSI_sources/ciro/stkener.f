      subroutine stkener(typener,enerx,ivalue)
      implicit real*8 (a-h,o-z)
      character*40 typener,typ
      dimension enerx(*)
      rewind 33
    1 read(33,end=3000)typ
      go to 1
3000  write(33)typener,(enerx(i),i=1,ivalue)
      rewind 33
      write(6,*)' *************************************'
      write(6,*) typener,' stockee sur la file 33'
      write(6,'(5e13.5)')(enerx(i),i=1,ivalue)
      write(6,*)' *************************************'
      return
      end
