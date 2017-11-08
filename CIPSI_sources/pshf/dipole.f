      subroutine dipole(scftyp)
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
c     common/hfden/da(5050),db(5050),ia(100)
      common/hfden/da(doas),db(doas),ia(doa)
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(dc),c(3,dc)
      common/iofile/ir,iw,ip,is,iq,ih,iv,if2,if3,if4,ioda(19)
c     common/dipxyz/xs(5050),ys(5050),zs(5050)
      common/dipxyz/xs(doas),ys(doas),zs(doas)
      data open /8huhf     /
      data fac /2.54158059d+00/
      if(ich.ne.0) return
c
c     ----- calculate dipole moment integrals -----
c
c
c     ----- electronic contribution to dipole moment -----
c
      dmx=-tracp(da,xs,num)
      dmy=-tracp(da,ys,num)
      dmz=-tracp(da,zs,num)
      if(scftyp.ne.open) go to 100
      dmx=-tracp(db,xs,num)+dmx
      dmy=-tracp(db,ys,num)+dmy
      dmz=-tracp(db,zs,num)+dmz
  100 continue
c
c     ----- nuclear contribution
c
      do 200 i=1,nat
      dmx=dmx+zan(i)*c(1,i)
      dmy=dmy+zan(i)*c(2,i)
      dmz=dmz+zan(i)*c(3,i)
  200 continue
      dipol=fac*sqrt(dmx*dmx+dmy*dmy+dmz*dmz)
      dmx=dmx*fac
      dmy=dmy*fac
      dmz=dmz*fac
      write(iw,9999) dipol,dmx,dmy,dmz
 9999 format(//,15x,'dipole moment',f10.4,' debye(s)',/,25x,'dmx',f10.4,
     1 /,25x,'dmy',f10.4,/,25x,'dmz',f10.4)
      return
      end
