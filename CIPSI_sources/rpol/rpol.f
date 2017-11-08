      program rcut
      include 'pshf.prm'
      character*26 timeday


      call fdate(timeday)
      write(6,8999) timeday
      call openf

      call set
      call efield
      call fdate(timeday)
      write(6,8999) timeday
 8999 format(/,80(1h*),/,8(10h   rpol   ),/, 10x,a25,/,
     * 80(1h*))
      stop
      end
