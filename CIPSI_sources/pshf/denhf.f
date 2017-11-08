      subroutine denhf(scftyp)
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      common/hfden/da(doas),db(doas),ia(doa)
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(dc),c(3,dc)
      common/iofile/ir,iw,ip,is,iq,ih,iv,if2,if3,if4,ioda(19)
      data open /8huhf     /
      do 50 i=1,num
   50 ia(i)=(i*(i-1))/2
      call reada(da,nx,6)
      if(scftyp.ne.open) return
      call reada(db,nx,9)
      return
      end
