      subroutine hfprop(scftyp)
      implicit real*8 (a-h,o-z)
      common/runopt/runtyp
      data uhf/8huhf     /
c
c   ----- calculate spin state for uhf wf.
c
      if(scftyp.eq.uhf) then
      call spin (sz,s2)
      end if
c
c     ----- dipole moment -----
c
      call denhf(scftyp)
      call dipole(scftyp)
c
c     ----- mulliken population analysis -----
c
      call denhf(scftyp)
      call mulken(scftyp)
c
c     ----- atomic spin density -----
c
      call denhf(scftyp)
      call spind(scftyp)
      return
      end
