      subroutine stktyp(typener,ispin,typsym,ivalue)
      implicit real*8 (a-h,o-z)
      character*40 typener,typ
      character  *3 typsym(*)
      dimension ispin(*)
      rewind 33
    1 read(33,end=3000)typ
      go to 1
3000  write(33)typener,(ispin(i),typsym(i),i=1,ivalue)
      rewind 33
      write(6,*)' *************************************'
      write(6,*) typener,' stockee sur la file 33'
      write(6,9876)(ispin(i),typsym(i),i=1,ivalue)
      write(6,*)' *************************************'
 9876 format(20(1x,i1,a3))
      return
      end
