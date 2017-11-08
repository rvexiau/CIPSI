C ***************************************************************
C calculate table of matrix elements for given l and radial function
C
C     <Y(l,m)f(r)|fpol(r)*fcut(r)|Y(l,m')f'(r)>
C
C   Y(l,m) : real,normalized spherical harmonics
C
C
C f(r) =  exp(-gama(r-ra)^2)
C f'(r)=  exp(-gamb(r-rb)^2)
C
C ysoft=TRUE ::  fcut = (1-exp(-(r/rcut)^2))^2 
C ysoft=FALSE::  fcut = theta(|r-rcut|)
C
C     value returned in mmtab(2l+1,2l+1,i)
C
C where fpol is one of the following according to ivec:
C
C ivec = 1: 1st integral (1/r^4 term) 
C ivec = 2: 2nd integral (x/r^3 term) 
C ivec = 3: 2nd integral (y/r^3 term) 
C ivec = 4: 2nd integral (z/r^3 term) 

C rp,ra,rb position of polarisable atom and of the two centers of gaussians
C this function uses the rtab table initialized by rint

      subroutine xyztrans(tab,rp,ra,rb,la,lb,gama,gamb,rcut,ysoft)
      implicit none

      real*8 tab(4,10,10)
      real*8 rp(3),ra(3),rb(3)
      integer la,lb
      real*8 gama,gamb
      real*8 rcut
      logical ysoft

      real*8 r4tab(28),r3tab(36)

      integer l1(3),l2(3)
      real*8 r1(3),r2(3),rsvec(3)
      real*8 rs
      real*8 gamma,coeff

      real*8 rmat(3,3)

      real*8 xyz0
 
      real*8 wigmat(2,7,7)

      integer i,j,l

      integer polyz(5), spherz(5)
      data polyz/0,1,4,10,20/, spherz/0,1,4,9,16/

      integer ltab(3,20)

C This table contains the angular momenta
C                                               
      data ltab
     &/ 0, 0, 0, 
     &  0, 0, 1, 
     &  1, 0, 0, 
     &  0, 1, 0,
     &  2, 0, 0, 
     &  0, 2, 0, 
     &  0, 0, 2, 
     &  1, 1, 0, 
     &  1, 0, 1, 
     &  0, 1, 1,
     &  3, 0, 0, 
     &  0, 3, 0, 
     &  0, 0, 3, 
     &  2, 1, 0, 
     &  2, 0, 1, 
     &  1, 2, 0, 
     &  0, 2, 1, 
     &  1, 0, 2, 
     &  0, 1, 2, 
     &  1, 1, 1/



C calculate center of gravity of the 2 gaussians and new exponent
      gamma=gama+gamb
      coeff=exp(-gama*gamb/gamma*
     &      ((ra(1)-rb(1))**2 +(ra(2)-rb(2))**2 +(ra(3)-rb(3))**2))
      do i=1,3
        rsvec(i) = (gama*ra(i)+gamb*rb(i))/gamma
      end do



C calculate rotation matrix, shift polarisable atom to origin
C and rotate rsvec to z - axis

      do i=1,3
        rsvec(i) = rsvec(i)-rp(i)
      end do
        
      call rot3d(la,lb,rsvec,rs,rmat,wigmat)
      do i=1,3
        r1(i)=0.0
        r2(i)=0.0
        do j=1,3
           r1(i) = r1(i) + rmat(i,j)*(ra(j)-rp(j))
           r2(i) = r2(i) + rmat(i,j)*(rb(j)-rp(j))
        end do
      end do

      l = max(la,lb)

C initialize integral tables
      call rint(l,rs,gamma,rcut,ysoft,r4tab,r3tab)

C fill table of matrix elements between simple cartesians
      do i=1,polyz(la+2)-polyz(la+1)
        l1(1) = ltab(1,i+polyz(la+1))
        l1(2) = ltab(2,i+polyz(la+1))
        l1(3) = ltab(3,i+polyz(la+1))
        do j=1,polyz(lb+2)-polyz(lb+1)
          l2(1) = ltab(1,j+polyz(lb+1))
          l2(2) = ltab(2,j+polyz(lb+1))
          l2(3) = ltab(3,j+polyz(lb+1))

          tab(1,i,j) =  coeff*xyz0(l1,l2,r1,r2,0,r4tab)
          tab(2,i,j) =  coeff*xyz0(l1,l2,r1,r2,1,r3tab)
          tab(3,i,j) =  coeff*xyz0(l1,l2,r1,r2,2,r3tab)
          tab(4,i,j) =  coeff*xyz0(l1,l2,r1,r2,3,r3tab)
        end do
      end do
C transform to real spherical harmonics and normalize the angular part
      call polytospher(la,lb,gama,gamb,tab)

C wigner rotation back
      call wigrot(la,lb,tab,wigmat,rmat)

      end

