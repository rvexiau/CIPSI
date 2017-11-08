      program moybig
      implicit real*8(a-h,o-x,z),logical*1(y)
      character*9 day,hour
      day='aujourdhui'
      hour='maintenant'
c      call date(day)
c      call time(hour)
c      write(6,9999)day,hour
 9999 format(/,80(1h*),/,8(10h  moyen   ),/, 10x,a10,5x,a10,/,
     * 80(1h*))
       call openf
       call moy
c      call date(day)
c      call time(hour)
      write(6,9998)day,hour
 9998 format(/,80(1h*),/,8(10h fin moyen),/, 10x,a10,5x,a10,/,
     * 80(1h*))
      stop
      end
