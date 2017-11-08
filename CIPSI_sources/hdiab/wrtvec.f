      subroutine wrtvec(v,nv,ndet)
      implicit real*8(a-h,o-z)
      dimension v(*)
      do  i=1,ndet
            WRITE(6,595) (v((j-1)*ndet+i),j=1,nv)
      enddo
595   format(10f12.6)
      return
      end
