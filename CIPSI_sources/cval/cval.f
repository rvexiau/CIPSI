      program cval
      implicit double precision (a-h,o-z)
      include 'pshf.prm'
      character *26 timeday
      logical ymono
      common/infnp/nprint
      common/iofile/ir,iw,ip,is,iq,ih,iv,if2,if3,if4,ioda(19)
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(dc),c(3,dc),
     1 ymono,mrec
      common/polar/calfa(dc)
      namelist /vpol/ calfa,zan,nprint,ymono,mrec
      if2=2
      if3=3
      ir=5
      iw=6
      is=8
      iq=9
      ih=10
      nprint=0
      call openf

      call fdate(timeday)
      write(6,8999) timeday
      call set
      read(ir,vpol)
      if (ymono) then
        write (6,*) 'pas de bielectroniques'
      end if
      do i=1,nat
      write(6,*)' atome ',i,' polarisabilite ',calfa(i),' charge nette',
     * zan(i)
      end do
      call polzer
      call polone
      if (.not.ymono) call poltwo
      call daclos

      call fdate(timeday)
      write(6,8999) timeday
8999  format(/,80(1h*),/,8(10h   cval   ),/, 10x,a25,/,80(1h*))
      stop
      end
