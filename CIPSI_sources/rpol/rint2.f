C ***************************************************************
C *  calculates integrals (ic=1) 
C I1 = Int dr dtheta  r^(l-2) Sin(th)**k
C      exp(-gamma(r^2 + rs^2-2r*rs*cos(th)) fcut(r,rcut)
C
C and (ic=2)
C
C I2 =  Int dr dtheta  r^(l-2) Sin(th)**k Cos(th)
C       exp(-gamma(r^2 + rs^2-2r*rs*cos(th)) fcut(r,rcut)
C
C ysoft=TRUE ::  fcut = (1-exp(-(r/rcut)^2))^2 
C ysoft=FALSE::  fcut = theta(|r-rcut|)
C

C this subroutine performs the r integration (Simpson rule)
C fills tables for l=0...lmax, k=1,l+1: r4tab(k l) and r3tab(k l)
C for 1/r**4 and r/r**3 integrals resp.
      subroutine rint(lmax,rs,gamma,rcut,ysoft,r4tab,r3tab)
      implicit none

      real*8 EPSI
      integer niter     
      parameter (EPSI=1.D-8,niter=50)

      integer lmax
      real*8 rs,gamma,rcut
      logical ysoft

      real*8 r4tab(28)
      real*8 r3tab(36)
      real*8 thetaint

      integer ictab4(28)
      integer ictab3(36)

      real*8 x
      integer i,k,l
      integer ix,nit
      real*8 r,rmin,rmax,rstep

      integer ic

C This tables contains a 1 or 2, if the corresponding value in rtab is needed
C in the calculation of the polarization energy. 1 means that integral of type
C I1 is stored, 2 of type I2, resp.; 0 in ictab: value will never be needed

      data ictab4 /1,2,0,1,0,1,2,0,2,0,
     &             1,0,1,0,1,2,0,2,0,2,
     &             0,1,0,1,0,1,0,1/

      data ictab3 /0,2,0,1,0,1,2,0,2,0,
     &             1,0,1,0,1,2,0,2,0,2,
     &             0,1,0,1,0,1,0,1,2,0,
     &             2,0,2,0,2,0/


      rmax = sqrt(-log(EPSI)/gamma)+rs
C estimate for integration boundary
      if (ysoft) then
         rmin = 0.0
      else
          rmin = rcut
      end if

C empty interval -- nothing to do
      if (rmin.ge.rmax) then
        do i=1,28
          r4tab(i)=0.0
        end do
        do i=1,36
          r3tab(i)=0.0
        end do
        return 
      end if

      rstep = (rmax-rmin)/(2.*(niter+1))
      if (int(rcut/rstep).lt.niter/4) then
C smaller step needed to increase precision at small r
            rstep = 4.*rcut/niter
	    nit = int((rmax-rmin)/(2.*rstep)) -1
            rstep = (rmax-rmin)/(2.*(nit+1))
      else
            nit = niter
      end if
      print *,'rint: # iterations: ',nit

C fill r4tab table

      ix = 1
      do l=0,2*lmax
        do k=1,l+1
          ic = ictab4(ix)
          if (ic.ne.0) then
c           print *,'rint: tab1',ix
            x = thetaint(ic,k,l,rmin,rs,gamma,rcut,ysoft)
     &        + thetaint(ic,k,l,rmin+rstep,rs,gamma,rcut,ysoft)*4.0
     &        + thetaint(ic,k,l,rmax,rs,gamma,rcut,ysoft)

            r = rmin+2.*rstep

            do i=1,nit
              x = x+thetaint(ic,k,l,r,rs,gamma,rcut,ysoft)*2.0
     &          + thetaint(ic,k,l,r+rstep,rs,gamma,rcut,ysoft)*4.0
              r = r + 2.*rstep
            end do
            r4tab(ix) = x*(rmax-rmin)/(6.*(nit+1))
          else
            r4tab(ix) = 9.999999d99
          end if
 
	  ix=ix+1
        end do
      end do

Cnow fill r3tab 

      ix = 1
      do l=0,2*lmax+1
        do k=1,l+1
          ic = ictab3(ix)
          if (ic.ne.0) then
c           print *,'rint: tab2',ix
            x = thetaint(ic,k,l+1,rmin,rs,gamma,rcut,ysoft)
     &        + thetaint(ic,k,l+1,rmin+rstep,rs,gamma,rcut,ysoft)*4.0
     &        + thetaint(ic,k,l+1,rmax,rs,gamma,rcut,ysoft)

            r = rmin+2.*rstep

            do i=1,nit
              x = x+thetaint(ic,k,l+1,r,rs,gamma,rcut,ysoft)*2.0
     &          + thetaint(ic,k,l+1,r+rstep,rs,gamma,rcut,ysoft)*4.0
              r = r + 2.*rstep
            end do
            r3tab(ix) = x*(rmax-rmin)/(6.*(nit+1))
          else
            r3tab(ix) = 9.999999d99
          end if
 
	  ix=ix+1
        end do
      end do

99    return
      end
