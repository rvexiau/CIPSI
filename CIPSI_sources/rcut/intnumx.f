      subroutine intnumx
c
c
c     ******************************************************************
c     ******************************************************************
c
c
c
      implicit integer(a-z)
      double precision vala,valai,valb,valbi,ga,gb
      parameter(u=20)
      common/nocas/casp
      common/iselec/ival,ivalf
      common/comvala/vala(0:u),valai(0:3,0:u)
      common/comvalb/valb(0:u),valbi(0:3,0:u)
      common/vintgab/ga,gb
c
c
c
 5000 format(1x,///,1x,'donnees : ',/,1x,9('*'),//)
 5001 format(1x,///,1x,72('*'),///)
 5003 format(1x,'calcul des integrales 1/r**4',/,1x,28('*'),//)
 5004 format(1x,'calcul des integrales x,y,z/r**3',/,1x,32('*'),//)
c
c
c
      call varcom
c
c
c
      if(casp.eq.0)then
                 if(ival.eq.1) call intgab
                 if(ivalf.eq.2) call intgabx
                   elseif(casp.eq.1)then
                        if(ival.eq.1) call intga(valb,valbi)
                        if(ivalf.eq.2) call intgax(valb,valbi)
                   elseif(casp.eq.2)then
                        if(ival.eq.1) call intga(vala,valai)
                        if(ivalf.eq.2) call intgax(vala,valai)
                   else
                                if(ival.eq.1) call caspgam
                                if(ivalf.eq.2) call caspgax
      endif
c
c
c
      return
      end
