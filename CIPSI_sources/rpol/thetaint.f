C calculates the integrand for rint

      real*8 function thetaint(ic,k,l,r,rs,gamma,rcut,ysoft)
      implicit none

      real*8 LN1R,EPSIO,EPSIE
c      parameter(LN1R=1.227947177,EPSIO=0.001,EPSIE=1.d-7)
      parameter(LN1R=1.227947177,EPSIO=0.001,EPSIE=1.d-8)

      integer ic,k,l
      real*8 r,rs,gamma,rcut
      logical ysoft

      real*8 bessi
      real*8 a,aa,em,ep
      integer i

      real*8 tab1(4,4),tab2(5,4)
      real*8 atab1(5,4),atab2(4,4)
      real*8 tablim1(4),tablim2(4)
      real*8 gamtab(4)

      data tab1
     &/   1.,   0.,   0.,  0.,  
     &   -2.,   2.,   0.,  0., 
     &   24., -24.,   8.,  0., 
     & -720., 720.,-288., 48./

      data tab2
     &/     1.,   -1.,    0.,   0.,   0.,  
     &     -6.,    6.,   -2.,   0.,   0.,
     &    120., -120.,   48.,  -8.,   0.,
     &  -5040., 5040.,-2160., 480., -48./
 
      data atab1
     &/         2.0, .3333333333, .1666666667, .0003968354, .55114D-5,
     &  1.333333333, .1333333333, .0047619048, .8818342D-4, .100108D-5,
     &  1.066666667, .0761904762, .0021164021, .3206670D-4, .308333D-6,
     &  .9142857143, .0507936508, .0011544012, .1480001D-4, .123333D-6/

      data atab2
     &/ .6666666667, .0666666667, .0023809524, .4409171D-4,
     &  .2666666667, .0190476192, .0005291005, .8016675D-5,
     &  .1523809524, .0084656085, .0001924002, .2466669D-5,
     &  .1015873016, .0046176046, .8880009D-4, .9866677D-6/

      data gamtab
C sqrt(pi) * 2**(k-1) * gamma(k-1/2)
     &/  3.14159265359,3.14159265359,9.42477796077,47.1238898038/

      data tablim1
     &/1.57079632679, 1.1780972451,  0.9817477042,  0.8590292412/
    
      data tablim2
     &/0.3926990817, 0.1963495408,  0.122718463, 0.08590292412/
    

      a = 2.*gamma*r*rs

      if (ic.eq.1) then
C ic = 1
        if ((k/2)*2.eq.k) then
C k even
          if (a**(k/2).lt.EPSIE) then
            thetaint = tablim1(k/2)*exp(-gamma*(r*r+rs*rs))
	  else
            thetaint = gamtab(k/2+1)*bessi(k/2,a)/a**(k/2) 
     &               *exp(-gamma*(r*r+rs*rs))
	  end if
        else
C k odd
          if (a**k.lt.EPSIO) then
	    aa = a*a
	    thetaint=atab1(1,k/2+1)+aa*(atab1(2,k/2+1) + 
     &        aa*(atab1(3,k/2+1)+aa*(atab1(4,k/2+1)+aa*atab1(5,k/2+1))))
            thetaint = thetaint*exp(-gamma*(r*r+rs*rs))
          else
            ep = tab1(k/2+1,k/2+1)
            em = -ep
            do i=k/2,1,-1
              em =  a*em - tab1(i,k/2+1) 
              ep = -a*ep + tab1(i,k/2+1)  
            end do
            thetaint = -(em*exp(-gamma*(r-rs)**2) + 
     &                  ep*exp(-gamma*(r+rs)**2))/a**k    
	  end if
        end if 
      else
C ic = 2
        if ((k/2)*2.eq.k) then
C k even
          if (a.lt.EPSIE) then
	    thetaint = a*tablim2(k/2)*exp(-gamma*(r*r+rs*rs))
	  else
	    thetaint = gamtab(k/2+1)*bessi(k/2+1,a)/a**(k/2)
     &              *exp(-gamma*(r*r+rs*rs))
	  end if
        else
C k odd
          if (a**k.lt.EPSIO) then
	    aa = a*a
	    thetaint = (atab2(1,k/2+1)+aa*(atab2(2,k/2+1)+
     &        aa*(atab2(3,k/2+1)+aa*(atab2(4,k/2+1) ))))*a
            thetaint = thetaint*exp(-gamma*(r*r+rs*rs))
	  else
            ep = tab2(k/2+2,k/2+1)
            em = -ep
            do i=k/2+1,1,-1
              em =  a*em - tab2(i,k/2+1) 
              ep = -a*ep + tab2(i,k/2+1)  
            end do
      
          thetaint = (em*exp(-gamma*(r-rs)**2) + 
     &                ep*exp(-gamma*(r+rs)**2))/a**(k+1)    
          end if 
	end if
      end if
      if (ysoft) then
        if (r/rcut.lt.EPSIE) then
          thetaint = thetaint*r**(l+2)/rcut**4
        else
          thetaint = thetaint*(1-exp(-(r/rcut)**2))**2*r**(l-2)
        end if
      else
        if (r/rcut.lt.EPSIE) then
          thetaint = 0
        else
          thetaint = thetaint*r**(l-2)
        end if
      end if

      return
      end
