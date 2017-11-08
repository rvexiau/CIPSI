C ***************************************************************
C Transformation from nonnormalized polynomial (overcomplete) 
C to normalized spherical harmonics basis
C

      subroutine polytospher(la,lb,gama,gamb,tab)
      implicit none

      integer la,lb
      real*8 gama,gamb
      real*8 tab(4,10,10)

      real*8 tac(10,7)
      integer iv,i
      real*8 c1,c2
      integer polyz(4), spherz(4)
      data polyz/1,3,6,10/, spherz/1,3,5,7/

      real*8 cp,cd1,cd2,cf1,cf2,cf3,cf4
      real*8 nos,nop,nod,nof

c nos = sqrt(1/8*sqrt(pi/2))*sqrt(4pi) 
      parameter(nos = 1.40310414554)

c nop = sqrt[3/32*sqrt(pi/2)]
      parameter(nop = 0.342780104985)

c nod = sqrt[15/128*sqrt(pi/2)]
      parameter(nod = 0.38323980804)

c nof = sqrt[105/512*sqrt(pi/2)] 
      parameter(nof = 0.506978612287)


cp = sqrt(3/(4pi))
      parameter(cp = 0.48860251190292)

cd1 = sqrt(5/(16pi))
      parameter(cd1 = 0.31539156525252)

cd2 = sqrt(15/(4pi))
      parameter(cd2 = 1.09254843059208)

cf1= sqrt(7/(16pi))
      parameter(cf1 = 0.37317633259012)

cf2= sqrt(21/(32pi))
      parameter(cf2 = 0.45704579946447)

cf3= sqrt(105/(16pi))
      parameter(cf3 = 1.44530572132028) 

cf4= sqrt(35/(32pi))
      parameter(cf4 = 0.59004358992664)

C normalization factor for gaussian functions
      if (la.eq.0) then
         c1 = nos/gama**0.75
      else if (la.eq.1) then
         c1 = nop/gama**1.25
      else if (la.eq.2) then
         c1 = nod/gama**1.75
      else if (la.eq.3) then
         c1 = nof/gama**2.25
      end if

      if (lb.eq.0) then
         c2 = nos/gamb**0.75
      else if (lb.eq.1) then
         c2 = nop/gamb**1.25
      else if (lb.eq.2) then
         c2 = nod/gamb**1.75
      else if (lb.eq.3) then
         c2 = nof/gamb**2.25
      end if
        
      do iv=1,4
        do i=1,polyz(la+1)
          if     (lb.eq.0) then 
             tac(i,1) = tab(iv,i,1)/c1
          else if(lb.eq.1) then
             tac(i,3) = cp*tab(iv,i,3)/c1
             tac(i,1) = cp*tab(iv,i,1)/c1
             tac(i,2) = cp*tab(iv,i,2)/c1
          else if(lb.eq.2) then
             tac(i,1) = cd1*(2*tab(iv,i,3)-tab(iv,i,1)-tab(iv,i,2))/c1
             tac(i,2) = cd2*tab(iv,i,5)/c1
             tac(i,3) = cd2*tab(iv,i,6)/c1
             tac(i,4) = cd2*(tab(iv,i,1)-tab(iv,i,2))/(2.*c1)
             tac(i,5) = cd2*tab(iv,i,4)/c1
          else if(lb.eq.3) then
             tac(i,1) = cf1*
     &                  (2*tab(iv,i,3)-3*tab(iv,i,5)-3*tab(iv,i,7))/c1
             tac(i,2) = cf2*(4*tab(iv,i,8)-tab(iv,i,1)-tab(iv,i,6))/c1
             tac(i,3) = cf2*(4*tab(iv,i,9)-tab(iv,i,2)-tab(iv,i,4))/c1
             tac(i,4) = cf3*(tab(iv,i,5)-tab(iv,i,7))/c1
             tac(i,5) = cf3*2*tab(iv,i,10)/c1
             tac(i,6) = cf4*(tab(iv,i,1)-3*tab(iv,i,6))/c1
             tac(i,7) = cf4*(tab(iv,i,2)-3*tab(iv,i,4))/c1
          end if 
        end do 

        do i=1,spherz(lb+1)
          if     (la.eq.0) then 
             tab(iv,1,i) = tac(1,i)/c2
          else if(la.eq.1) then
             tab(iv,3,i) = cp*tac(3,i)/c2
             tab(iv,1,i) = cp*tac(1,i)/c2
             tab(iv,2,i) = cp*tac(2,i)/c2
          else if(la.eq.2) then
             tab(iv,1,i) = cd1*(2*tac(3,i)-tac(1,i)-tac(2,i))/c2
             tab(iv,2,i) = cd2*tac(5,i)/c2
             tab(iv,3,i) = cd2*tac(6,i)/c2
             tab(iv,4,i) = cd2*(tac(1,i)-tac(2,i))/(2.*c2)
             tab(iv,5,i) = cd2*tac(4,i)/c2
          else if(la.eq.3) then
             tab(iv,1,i) = cf1*(2*tac(3,i)-3*tac(5,i)-3*tac(7,i))/c2
             tab(iv,2,i) = cf2*(4*tac(8,i)-tac(1,i)-tac(6,i))/c2
             tab(iv,3,i) = cf2*(4*tac(9,i)-tac(2,i)-tac(4,i))/c2
             tab(iv,4,i) = cf3*(tac(5,i)-tac(7,i))/c2
             tab(iv,5,i) = cf3*2*tac(10,i)/c2
             tab(iv,6,i) = cf4*(tac(1,i)-3*tac(6,i))/c2
             tab(iv,7,i) = cf4*(tac(2,i)-3*tac(4,i))/c2
          end if 
        end do 

      end do

      
      end

