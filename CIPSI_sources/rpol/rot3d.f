C ***************************************************************
C
C calculates rotation matrix
C for 3d-rotation of real space vectors        -> rmat
C 
C                                         (r(1))   ( 0 )
C  D = D_phiD_theta with i the property D (r(2)) = ( 0 )
C                                         (r(3))   (|r|)
C
C and the wigner rotation matrices for given la,lb -> wigmat
C
C form of wigmat:
C
C     ( D00  D01+(-1)^1D0,-1  D01-(-1)^1D0,-1 D02+(-1)^2D02  D02-(-1)^2D02,-2  
C     ( D10  D11+(-1)^1D1,-1  D11-(-1)^1D1,-1 D12+(-1)^2D12  D12-(-1)^2D12,-2  
C D = ( D20  D21+(-1)^1D2,-1  D21-(-1)^1D2,-1 D22+(-1)^2D22  D22-(-1)^2D22,-2  
C     ( D30  D31+(-1)^1D3,-1  D31-(-1)^1D3,-1 D32+(-1)^2D32  D32-(-1)^2D32,-2  
C     ( D40  D41+(-1)^1D4,-1  D41-(-1)^1D4,-1 D42+(-1)^2D42  D42-(-1)^2D42,-2  
C 
C
C returns also length of r vector in rs

      subroutine rot3d(la,lb,r,rs,rmat,wigmat)
      implicit none

      integer la,lb
      real*8 r(3) 
      real*8 rmat(3,3)
      real*8 rs
      real*8 wigmat(2,7,7)


      integer il,i,j,l
      real*8 rab
      real*8 mat(7,7)
      real*8 sini,cosi,cosold,sinphi,cosphi
      real*8 sign

      real*8 costh,cos2th,costh2
      real*8 sinth,sin2th,sinth2
      real*8 rsinth2,rcosth2,cos3th2,sin3th2 

      real*8 EPSI
c      parameter(EPSI=1.D-6)
      parameter(EPSI=1.D-8)
      real*8 SQRT2
      parameter(SQRT2=1.41421356237)
    
      rab = sqrt(r(1)*r(1)+r(2)*r(2))
      rs = sqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))

      if (rab.lt.EPSI) then
C no rotation necessary
        cosphi = 1.0
        sinphi = 0.0

C rotation matrix is unity
        rmat(1,2) = 0.
        rmat(1,3) = 0.
        rmat(2,1) = 0.
        rmat(2,2) = 1.
        rmat(2,3) = 0.
        rmat(3,1) = 0.
        rmat(3,2) = 0.
        if (r(3).ge.0.0) then
          costh = 1.
          rmat(1,1) = 1. 
          rmat(3,3) = 1.
        else
          costh  = -1.
          rmat(1,1) = -1. 
          rmat(3,3) = -1.
        end if
      else
        cosphi = r(1)/rab
        sinphi = r(2)/rab
        costh  = r(3)/rs

        rmat(1,1) =  r(1)*r(3)/(rab*rs)
        rmat(1,2) =  r(2)*r(3)/(rab*rs)
        rmat(1,3) = -rab/rs
        rmat(2,1) = -r(2)/rab
        rmat(2,2) =  r(1)/rab
        rmat(2,3) =  0.
        rmat(3,1) =  r(1)/rs
        rmat(3,2) =  r(2)/rs
        rmat(3,3) =  r(3)/rs
      end if

C now set up wigner matrix for complex spherical harmonics

      l = la
c     costh=1.0
      do il=1,2
************************************
*               l=0                *
************************************
        if (l.eq.0) then       
          mat(1,1) = 1.0
************************************
*               l=1                *
************************************
        else if (l.eq.1) then
          sinth = -sqrt(1.0-costh*costh)

          mat(1,1) =  (1.0+costh)/2.0
          mat(1,2) =  sinth/sqrt(2.0)
          mat(1,3) =  (1.0-costh)/2.0
          mat(2,1) = -mat(1,2)
          mat(2,2) =  costh
          mat(2,3) =  mat(1,2)
          mat(3,1) =  mat(1,3)
          mat(3,2) = -mat(1,2)
          mat(3,3) =  mat(1,1)
 
************************************
*               l=2                *
************************************
        else if (l.eq.2)  then

          sinth = -sqrt(1.0-costh*costh)
C costh2 = cos(th/2)**2
          costh2 = (1.0+costh)/2.0
C sinth2 = sin(th/2)**2
          sinth2 = (1.0-costh)/2.0
C cos2th = cos(2 th)    
          cos2th = 2.0*costh*costh-1.0
C sin2th = sin(2 th)    
          sin2th = 2.0*sinth*costh


          mat(1,1) =  costh2*costh2
          mat(1,2) =  costh2*sinth 
          mat(1,3) =  sqrt(3./8.)*sinth*sinth 
          mat(1,4) =  sinth2*sinth 
          mat(1,5) =  sinth2*sinth2 
  
          mat(2,1) = -mat(1,2)
          mat(2,2) =  (costh*(2*costh+1)-1.)/2.0 
          mat(2,3) =  sqrt(3./8.)*sin2th 
          mat(2,4) =  (costh-cos2th)/2.0 
          mat(2,5) =  mat(1,4) 

          mat(3,1) =  mat(1,3) 
          mat(3,2) = -mat(2,3) 
          mat(3,3) =  (1. +3.*cos2th)/4.0
          mat(3,4) =  mat(2,3) 
          mat(3,5) =  mat(1,3) 

          mat(4,1) = -mat(1,4) 
          mat(4,2) =  mat(2,4) 
          mat(4,3) = -mat(2,3)
          mat(4,4) =  mat(2,2) 
          mat(4,5) =  mat(1,2) 

          mat(5,1) =  mat(1,5) 
          mat(5,2) = -mat(1,4) 
          mat(5,3) =  mat(1,3)
          mat(5,4) = -mat(1,2) 
          mat(5,5) =  mat(1,1) 

************************************
*               l=3                *
************************************
        else if (l.eq.3) then

          sinth = -sqrt(1.0-costh*costh)

C costh2 = cos(th/2)**2
          costh2 = (1.0+costh)/2.0

C sinth2 = sin(th/2)**2
          sinth2 = (1.0-costh)/2.0

C cos2th = cos(2 th)    
          cos2th = 2.0*costh*costh-1.0

C sin2th = sin(2 th)    
          sin2th = 2.0*sinth*costh

c rsinth2 = sin(th/2)
          rsinth2 = -sqrt(sinth2)

c rcosth2 = cos(th/2)
          rcosth2 = sqrt(costh2)

c cos3th2 = cos(3*th/2)
          cos3th2 = costh*rcosth2-sinth*rsinth2

c sin3th2 = sin(3*th/2)
          sin3th2 = costh*rsinth2+sinth*rcosth2

          mat(1,1) =  costh2*costh2*costh2
          mat(1,2) =  sqrt(1.5)*costh2*costh2*sinth
          mat(1,3) =  sqrt(15.)*costh2*costh2*sinth2 
          mat(1,4) =  sqrt(5.)/4.*sinth*sinth*sinth
          mat(1,5) =  sqrt(15.)*costh2*sinth2*sinth2 
          mat(1,6) =  sqrt(1.5)*sinth2*sinth2*sinth
          mat(1,7) =  sinth2*sinth2*sinth2
 
          mat(2,1) = -mat(1,2)
          mat(2,2) =  costh2*costh2*(-2.+3.*costh)
          mat(2,3) =  sqrt(5./8)*rcosth2*costh2*
     &                       (3*sin3th2-5*rsinth2) 
          mat(2,4) =  sqrt(15./8)*costh*sinth*sinth
          mat(2,5) =  sqrt(5./8)*rsinth2*sinth2*
     &                       (3*cos3th2+5*rcosth2) 
          mat(2,6) =  sinth2*sinth2*(2.+3.*costh)
          mat(2,7) =  mat(1,6)

          mat(3,1) =  mat(1,3)
          mat(3,2) = -mat(2,3)
          mat(3,3) =  costh2*(costh*(15.*costh-10.)-1.)/4.
          mat(3,4) =  sqrt(3.)/4.*sinth*(5*costh*costh-1.)
          mat(3,5) =  sinth2*(costh*(15.*costh+10.)-1.)/4.
          mat(3,6) =  mat(2,5)
          mat(3,7) =  mat(1,5)

          mat(4,1) = -mat(1,4)   
          mat(4,2) =  mat(2,4)
          mat(4,3) = -mat(3,4)
          mat(4,4) =  costh*(2.5*costh*costh-1.5)
          mat(4,5) = -mat(4,3)
          mat(4,6) =  mat(4,2) 
          mat(4,7) = -mat(4,1)

          mat(5,1) =  mat(1,5)
          mat(5,2) = -mat(2,5) 
          mat(5,3) =  mat(3,5)
          mat(5,4) = -mat(3,4)
          mat(5,5) =  mat(3,3)
          mat(5,6) =  mat(2,3)
          mat(5,7) =  mat(3,1)

          mat(6,1) = -mat(1,6)
          mat(6,2) =  mat(2,6) 
          mat(6,3) = -mat(3,6)
          mat(6,4) =  mat(2,4)
          mat(6,5) = -mat(2,3)
          mat(6,6) =  mat(2,2)
          mat(6,7) =  mat(1,2)

          mat(7,1) =  mat(1,7)
          mat(7,2) = -mat(1,6) 
          mat(7,3) =  mat(1,5)
          mat(7,4) = -mat(1,4)
          mat(7,5) =  mat(1,3)
          mat(7,6) = -mat(1,2)
          mat(7,7) =  mat(1,1)
        else
          write (6,*) 'rotmat: l>3 not supported, l=',l
          stop
        end if

C now construct transformation matrix for real spherical harmonics
        sign = -1
        wigmat(il,1,1) = mat(l+1,l+1)
        do j=1,l
          wigmat(il,1,2*j)   =(mat(l+1,l+1+j)+sign*mat(l+1,l+1-j))/SQRT2
          wigmat(il,1,2*j+1) =0.0
          sign = -sign
        end do

        sini = -sinphi
        cosi =  cosphi

        do i=1,l
          sign = -1
          wigmat(il,2*i  ,1) = cosi*mat(l+1+i,l+1)*SQRT2
          wigmat(il,2*i+1,1) =-sini*mat(l+1+i,l+1)*SQRT2
          do j=1,l
             wigmat(il,2*i  ,2*j  ) = 
     &               cosi*(mat(l+1+i,l+1+j)+sign*mat(l+1+i,l+1-j))
             wigmat(il,2*i  ,2*j+1) = 
     &               sini*(mat(l+1+i,l+1+j)-sign*mat(l+1+i,l+1-j))
             wigmat(il,2*i+1,2*j  ) = 
     &              -sini*(mat(l+1+i,l+1+j)+sign*mat(l+1+i,l+1-j))
             wigmat(il,2*i+1,2*j+1) = 
     &               cosi*(mat(l+1+i,l+1+j)-sign*mat(l+1+i,l+1-j))
             sign = -sign
          end do
C cos(m phi) -> cos((m+1) phi)
          cosold=cosi
          cosi = cosi*cosphi+sini*sinphi
          sini = sini*cosphi-cosold*sinphi

        end do
        l = lb
      end do    
      return
      end





