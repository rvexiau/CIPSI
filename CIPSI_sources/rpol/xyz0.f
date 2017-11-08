C ***************************************************************
C *  calculates integrals 
C   Int dr dphi dtheta r^2 * sin(th)* (x-x1)^lx1 * (y-y1)^ly1 * (z-z1)^lz1 * 
C                                     (x-x2)^lx2 * (y-y2)^ly2 * (z-z2)^lz2 * 
C      r^(l-2)*exp(-gamma(r^2 + rs^2-2r*rs*cos(th)) fcut(r,rcut)
C
C   and
C
C   Int dr dphi dtheta r^2 * sin(th)* (x-x1)^lx1 * (y-y1)^ly1 * (z-z1)^lz1 * 
C                 (x y z)             (x-x2)^lx2 * (y-y2)^ly2 * (z-z2)^lz2 * 
C      *r^(l-1)*exp(-gamma(r^2 + rs^2-2r*rs*cos(th)) fcut(r,rcut)
C
C      with x = r cos(ph) sin(th); y = r sin(ph) sin(th); z = r cos(th)
C
C ysoft=TRUE ::  fcut = (1-exp(-(r/rcut)^2))^2 
C ysoft=FALSE::  fcut = theta(|r-rcut|)
C
C ivec = 0: 1st integral (1/r^4 term)  use r4tab !
C ivec = 1: 2nd integral (x/r^3 term)  use r3tab!
C ivec = 2: 2nd integral (y/r^3 term)      "
C ivec = 3: 2nd integral (z/r^3 term)      "

C this function uses the rtab table initialized by rint

      real*8 function xyz0(l1,l2,r1,r2,ivec,rtab)
      implicit none

      integer l1(3),l2(3)
      real*8 r1(3),r2(3)
      integer ivec
      real*8 rtab(36)

      integer i1x,i1y,i1z
      integer i2x,i2y,i2z
      integer ix,iy,iz
      real*8 f1x,f1y,f1z
      real*8 f2x,f2y,f2z

      real*8 xyzint

      real*8 binom(7,7)

C This table contains the binomial coefficients (i)
C                                               (j)
      data binom
     &/ 1., 0., 0., 0., 0., 0., 0., 
     &  1., 1., 0., 0., 0., 0., 0.,
     &  1., 2., 1., 0., 0., 0., 0.,
     &  1., 3., 3., 1., 0., 0., 0.,
     &  1., 4., 6., 4., 1., 0., 0.,
     &  1., 5.,10.,10., 5., 1., 0.,
     &  1., 6.,15.,20.,15., 6., 1./ 

      if      (ivec.eq.0) then
       ix = 0
       iy = 0
       iz = 0
      else if (ivec.eq.1) then
       ix = 1
       iy = 0
       iz = 0
      else if (ivec.eq.2) then
       ix = 0
       iy = 1
       iz = 0
      else if (ivec.eq.3) then
       ix = 0
       iy = 0
       iz = 1
      else 
        write(6,*) 'Error xyz0: invalid ivec: ',ivec
        stop
      end if

      xyz0 = 0.0
C lx and ly must be even, else xyint is zero
      f1x = 1.
      do i1x=l1(1),0,-1
        f1y = 1.
        do i1y=l1(2),0,-1
          f1z = 1.
          do i1z=l1(3),0,-1
            f2x = 1.
            do i2x=l2(1),0,-1
              f2y = 1.
              do i2y=l2(2),0,-1
                f2z = 1.
                do i2z=l2(3),0,-1

                  xyz0 = xyz0 +
     &              xyzint(i1x+i2x+ix,i1y+i2y+iy,i1z+i2z+iz,rtab)*
     &              binom(i1x+1,l1(1)+1)*binom(i1y+1,l1(2)+1)*
     &              binom(i1z+1,l1(3)+1)*f1x*f1y*f1z*
     &              binom(i2x+1,l2(1)+1)*binom(i2y+1,l2(2)+1)*
     &              binom(i2z+1,l2(3)+1)*f2x*f2y*f2z

                  f2z = -f2z*r2(3)
                end do
                f2y = -f2y*r2(2)
              end do
              f2x = -f2x*r2(1)
            end do
            f1z = -f1z*r1(3)
          end do
          f1y = -f1y*r1(2)     
        end do
        f1x = -f1x*r1(1)
      end do

      end

