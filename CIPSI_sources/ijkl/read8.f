      subroutine read8(xx,ix,lsize)
      implicit real*8 (a-h,o-x,z)
      dimension xx(lsize),ix(lsize)
      read(8) xx ,ix
      return
      end
