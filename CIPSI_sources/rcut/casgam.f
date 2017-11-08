      subroutine casgam(a)
c
c
c     cette subroutine calcule la fonction gamma entre les bornes "a" et
c     "l'infini". cette integrale est obtenue en soustrayant l'integrale
c     entre les bornes "zero" et "a" calculee par les points de gauss de
c     l'integrale entre les bornes "zero" et "l'infini" que l'on peut ob-
c     tenir dans les tables (a l'aide de formules mathematiques simples).
c
c
c     ..................................................................
c
c
c
      implicit  double precision(a-h,o-z)
      double precision intainf
      integer ent,sor
      common/numeriq/xg(64),hg(64)
      common/valpi/pi
      common/valeurz/z
      common/puissan/puist
      common/integr/intainf
      common/comvar/cte(2),pt
      common/lecrit/ent,sor
      data np2/64/
c
c
c     calcul de l'integrale de la fonction gamma de "0" a "l'infini".
c     ---------------------------------------------------------------
c
c
c     on calcule gamma(z), il faut alors determiner si z est un entier
c     ou un nombre fractionnaire afin d'appliquer les formules :
c          gamma(n+1)=1*2*3...(n-1)n=n! ou
c          gamma(n+1/2)=(1*3*5*7...(2n-1)/2**n)*gamma(1/2)
c
c     gam=valeur numerique de la fonction gamma(n+1) ou gamma(n+1/2).
c
c     z est utilise dans gamma(z).
c
c
c
      if(z.lt.1.0d0) go to 10
      val=z/(dint(z))
      dval=dabs(val-1.0d0)
      if(dval.gt.1.0d-8) go to 10
      go to 75
c
c
c     cas ou z est un nombre fractionnaire :
c
 10   gam1s2=dsqrt(pi)
      vn=z-0.5d0
      gamfrac=1.0d0
      do 50 i=1,(2*idnint(vn)-1),2
      gamfrac=gamfrac*dfloat(i)
 50   continue
      gamfrac=(gamfrac*gam1s2)/2.0d0**vn
      gam=gamfrac
      go to 150
c
c
c     cas ou z est un nombre entier :
c
 75   m=z-1.0d0
      gament=1.0d0
      do 100 j=1,m
      gament=gament*dfloat(j)
 100  continue
      gam=gament
c
c
c
c     calcul de l'integrale entre les bornes "0" et "a" par la methode
c     des points de gauss.
c     ------------------------------------------------------------------
c
c
c     vint=valeur numerique de la fonction gamma evaluee entre les bornes
c          0 et a.
c
c
 150  vint=0.0d0
      if(dabs(a).lt.1.0d-8) go to 250
      bma=a*0.5d0
      bpa=bma
      do 200 i=1,np2
      x=bma*xg(i)+bpa
      t=(x**puist)*dexp(-x*x)
      t=t*hg(i)
      vint=vint+t
 200  continue
      vint=vint*bma*2.0d0
c
c
c
 250  intainf=gam-vint
c
c
c
      return
      end
