      subroutine intab(a)
c
c
c     cette subroutine calcule l'integrale :
c               2*int(a,b)[u**(2z-1)*exp(-u**2) du ] entre les bornes
c     a et b par la methode de gauss-legendre.
c
c     la valeur numerique de l'integrale=intainf.
c
c     .................................................................
c
c
      implicit  double precision (a-h,o-z)
      double precision intainf
      integer ent,sor
      common/numeriq/xg(64),hg(64)
      common/borneb/b
      common/puissan/puist
      common/integr/intainf
      common/comvar/cte(2),pt
      common/lecrit/ent,sor
      data np2/64/
c
c
c
      bma=0.5d0*(b-a)
      bpa=0.5d0*(b+a)
c
c
c     vintab : integrale entre les bornes a et b.
c     ..........................................
c
c
      vintab=0.0d0
      do 100 i=1,np2
      x=bma*xg(i)+bpa
      t=(x**puist)*dexp(-x*x)
      t=t*hg(i)
      vintab=vintab+t
 100  continue
c
c
      vintab=2.0d0*vintab*bma
      intainf=vintab
c
c
      return
      end
