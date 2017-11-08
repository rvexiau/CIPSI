      subroutine rdener(typener,enerx,ivalue)
      implicit real*8 (a-h,o-z)
      character*40 typener,typ
      dimension enerx(*),ener(1000)
      rewind 33
    1 read(33,end=3000)typ,jvalue,(ener(i),i=1,jvalue)
      if(typ.ne.typener)go to 1
      do i=1,ivalue
	enerx(i)=ener(i)
      end do
      return
3000  write(6,*)typener,'inexistant sur la file 33' 
c      stop
      end
