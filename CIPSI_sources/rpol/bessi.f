Calculates I_0(x) source: Abramovitz, Segun p. 378
      real*8 function bessi0(x)
      implicit none
      real*8 x
      real*8 t

      integer i
      real*8 ax

      real*8 tab1(7),tab2(9)

      data tab1
     &/1.0, 3.5156229, 3.0899424, 1.2067492, 0.2659732,
     & 0.0360768, 0.0045813/ 
      data tab2
     &/0.39894228, 0.01328592, 0.00225319, -.00157565, 0.00916281,
     & -.02057706, 0.02635537, -.01647633, 0.00392377/
      
      ax = abs(x)
      t=ax/3.75
      if(t.lt.1.0) then
        t = t*t
        bessi0 = t*tab1(7)
        do i=6,1,-1
	  bessi0= bessi0*t+tab1(i)
        end do
      else
	t=1./t
	bessi0 = t*tab2(9)
	do i=8,1,-1
	  bessi0= bessi0*t+tab2(i)
        end do
	bessi0 = bessi0*dexp(ax)/dsqrt(ax)
      end if
      return
      end

Calculates I_1(x) source: Abramovitz, Segun p. 378
      real*8 function bessi1(x)
      implicit none
      real*8 x
      real*8 t

      integer i
      real*8 ax
      real*8 tab1(7),tab2(9)
      data tab1
     &/0.5, 0.87890594, 0.51498869, 0.15084934, 0.02658733,
     & 0.00301532, 0.00032411/ 
      data tab2
     &/0.39894228, -.03988024, -.00362018, 0.00163801, -.01031555,
     & 0.02282967, -.02895312, 0.01787654, -.00420059/
      
      ax = abs(x)
      t=ax/3.75
      if(t.lt.1.0) then
        t = t*t
        bessi1 = t*tab1(7)
        do i=6,1,-1
	  bessi1= bessi1*t+tab1(i)
        end do
	bessi1 = bessi1*x
      else
	t=1./t
	bessi1 = t*tab2(9)
	do i=8,1,-1
	  bessi1= bessi1*t+tab2(i)
        end do
	bessi1 = bessi1*exp(ax)/sqrt(ax)
	if (x.lt.0) bessi1=-bessi1
      end if
      return
      end

      real*8 function bessi(n,x)
      integer n
      real*8 x
      real*8 bessi0,bessi1

      real*8 tab2(4),tab3(4),tab4(4)
      data tab2/.125,.0104166667,0.0003255208,0.54253d-5/
      data tab3/.0208333333,.0013020833,.325521d-4,.4521d-6/
      data tab4/.0026041667,.0001302083,.27127d-5,.3229d-8/

      if (n.eq.0) then
	bessi= bessi0(x)
      else if (n.eq.1) then
	bessi= bessi1(x)
      else if (n.eq.2) then
	if (abs(x).lt.0.5) then
	  xx = x*x
	  bessi=(((xx*tab2(4)+tab2(3))*xx+tab2(2))*xx+tab2(1))*xx 
	else
	  bessi = -2./x*bessi1(x)+bessi0(x)
	end if

      else if (n.eq.3) then
	xx = x*x
	if (abs(x).lt.0.7) then
	  bessi=(((xx*tab3(4)+tab3(3))*xx+tab3(2))*xx+tab3(1))*xx*x 
	else
	  bessi = (8./xx+1.)*bessi1(x)-4./x*bessi0(x)
	end if

      else if (n.eq.4) then
	xx = x*x
	if (abs(x).lt.1.0) then
	  bessi=(((xx*tab4(4)+tab4(3))*xx+tab4(2))*xx+tab4(1))*xx*xx
	else
	  bessi = -(48/xx+8)/x*bessi1(x)+(1.+24./xx)*bessi0(x)
	end if

      else 
	write(6,*) 'Error bessi: order',n,' greater then 4'
	stop
      end if
      return
      end 

c      program tesbes
c
c      implicit none
c     real*8 bessi
c      real*8 x
c      integer n
c     
c      read (*,*) n,x
c      print *,' I(',n,' , ',x,' ) =  ',bessi(n,x)
c      
c      end
