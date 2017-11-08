C ***************************************************************
C *  calculates integrals 
C   Int dr dphi dtheta r^2 * sin(th)* x^lx * y^ly * z^lz * 
C      exp(-gamma(r^2 + rs^2-2r*rs*cos(th)) fcut(r,rcut)
C
C      with x = r cos(ph) sin(th); y = r sin(ph) sin(th); z = r cos(th)
C
C ysoft=TRUE ::  fcut = (1-exp(-(r/rcut)^2))^2 
C ysoft=FALSE::  fcut = theta(|r-rcut|)
C

C this function uses the rtab table initialized by rint

      real*8 function xyzint(lx,ly,lz,rtab)
      implicit none

      integer lx,ly,lz
      real*8 rtab(28)

      real*8 PI512 
      parameter (PI512=0.613592315154d-2)

      real*8 ph
      integer phitab(8,8)

C This table contains 512/Pi*Int_0^2Pi dphi Sin^k(phi) Cos^l(phi)

      data phitab
     &/ 1024,   0, 512,   0, 384,   0, 320,  0,
     &     0,   0,   0,   0,   0,   0,   0,  0,
     &   512,   0, 128,   0,  64,   0,  40,  0,
     &     0,   0,   0,   0,   0,   0,   0,  0,
     &   384,   0,  64,   0,  24,   0,  12,  0,
     &     0,   0,   0,   0,   0,   0,   0,  0,
     &   320,   0,  40,   0,  12,   0,   5,  0,
     &     0,   0,   0,   0,   0,   0,   0,  0/

      xyzint = 0
      if (phitab(lx+1,ly+1).eq.0) then
        return
      else
        ph = phitab(lx+1,ly+1)*PI512
      end if

      if      (lz.eq.0) then
C k = lx+ly+1; l= lx+ly
        xyzint = ph*rtab((lx+ly+1)*(lx+ly+2)/2)

      else if (lz.eq.1) then
C k = lx+ly+1; l= lx+ly+1
        xyzint = ph*rtab((lx+ly+1)*(lx+ly+4)/2)

      else if (lz.eq.2) then
C k = lx+ly+1; l= lx+ly+2 and k=lx+ly+3; l=lx+ly+2
        xyzint = ph*(rtab((lx+ly+2)*(lx+ly+5)/2 - 1) - 
     &               rtab((lx+ly+2)*(lx+ly+5)/2 + 1))

      else if (lz.eq.3) then
C k = lx+ly+1; l= lx+ly+3 and k=lx+ly+3; l=lx+ly+3
        xyzint = ph*(rtab((lx+ly+3)*(lx+ly+6)/2 - 2) - 
     &               rtab((lx+ly+3)*(lx+ly+6)/2    ))

      else if (lz.eq.4) then
C k = lx+ly+1; l= lx+ly+4, k=lx+ly+3; l=lx+ly+4 and k=lx+ly+5; l=lx+ly+4
        xyzint = ph*(rtab((lx+ly+3)*(lx+ly+8)/2 - 1) - 
     &            2.*rtab((lx+ly+3)*(lx+ly+8)/2 + 1) + 
     &               rtab((lx+ly+3)*(lx+ly+8)/2 + 3))

      else if (lz.eq.5) then
C k = lx+ly+1; l= lx+ly+5, k=lx+ly+3; l=lx+ly+5 and k=lx+ly+5; l=lx+ly+5
        xyzint = ph*(rtab((lx+ly+4)*(lx+ly+9)/2 - 2) - 
     &            2.*rtab((lx+ly+4)*(lx+ly+9)/2    ) + 
     &               rtab((lx+ly+4)*(lx+ly+9)/2 + 2))

      else if (lz.eq.6) then
C k = lx+ly+1; l= lx+ly+6, k=lx+ly+3, k=lx+ly+5 and k=lx+ly+7
        xyzint = ph*(rtab((lx+ly+4)*(lx+ly+11)/2    ) - 
     &            3.*rtab((lx+ly+4)*(lx+ly+11)/2 + 2) + 
     &            3.*rtab((lx+ly+4)*(lx+ly+11)/2 + 4) - 
     &               rtab((lx+ly+4)*(lx+ly+11)/2 + 6))

      else if (lz.eq.7) then
C k = lx+ly+1; l= lx+ly+7, k=lx+ly+3, k=lx+ly+5 and k=lx+ly+7
        xyzint = ph*(rtab((lx+ly+5)*(lx+ly+12)/2 - 1) - 
     &            3.*rtab((lx+ly+5)*(lx+ly+12)/2 + 1) + 
     &            3.*rtab((lx+ly+5)*(lx+ly+12)/2 + 3) - 
     &               rtab((lx+ly+5)*(lx+ly+12)/2 + 5))
      end if

      end

