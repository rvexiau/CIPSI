      subroutine chelec
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      character*40 typener
      common/iofile/ir,iw,ip,is,iq,ih,iv
      common/times/ti,tx,tim
      common/infoa/nat,ich,mul,numscf,nx,ne,na,nb,zan(dc),c(3,dc),en
      common/dipxyz/xs(doas),ys(doas),zs(doas)
      common/elchp/ elec(3)
      common/output/nprint,itol,icut,normf,normp,minp
      dimension h(doas)
c
c    modification des integrales monoelectroniques dans un champ electri
c    uniforme
      elec2=dsqrt(elec(1)**2+elec(2)**2+elec(3)**2)
      if(elec2.gt.1.d-12) then
        call reada (h,nx,3)
        do i=1,nx
        h(i)=h(i)+elec(1)*xs(i)
        h(i)=h(i)+elec(2)*ys(i)
        h(i)=h(i)+elec(3)*zs(i)
        end do
        call wrtda(h,nx,3)
 

c   contribution nucleaire dans le champ electrique <-- taken into account in mole.f
c        do i=1,nat
c        rr=elec(1)*c(1,i)+elec(2)*c(2,i)+elec(3)*c(3,i)
c        en=en-zan(i)*rr
c        end do
      endif
      return
      end
