      subroutine intinf
c
c
c     cette subroutine calcule l'integrale de a a l'infini dans le cas
c     d'une puissance de u < 0 ou = 0 lors de l'evaluation de l'inte-
c     grale :
c
c                 2*int(a,infini)[u**(2z-1)*exp(-u**2) du]
c
c
c     ..................................................................
c
c
c
      implicit double precision (a-h,o-z)
      integer ent,sor
      parameter (npar=50)
      common/numeriq/xg(64),hg(64)
      common/borneb/b
      common/valpi/pi
      common/puissan/puist
      common/lecrit/ent,sor
      data np2/64/
c
c
 5000 format(1x,f25.15)
c
c
c
c     il s'agit de determiner 2*int(a,infini)[u**(2z-1)*exp(-u**2) du]
c     de la facon suivante :
c
c     2*int(a,infini)[u**(2z-1)*exp(-u**2) du]=2*int(a,b)+2*int(b,inf)
c
c
c
c     determination de la borne b:
c     ...........................
c
c     2*int(b,inf)[exp(-u**2)/u**(k) du <
c                          (2/b**(k))*int(b,infini)[exp(-u**2)du],
c                             avec k=(2z-1)
c     pour differentes valeurs de b, on evalue :
c     int(b,infini)=int(0,inf)[exp(-u**2)du] - int(0,b)[exp(-u**2)du],
c     la premiere integrale est egale a 1/2*(sqrt(pi)) et la seconde est
c     evaluee par la methode des points de gauss.
c     on obtient la valeur de b quand
c           (2/b**(k))*int(b,infini)[exp(-u**2)du] < 10**(-10).
c
c
c
c     zerinf : integrale de 0 a l'infini.
c     ..................................
c
      zerinf=0.5d0*dsqrt(pi)
c
c
c
      do 500 j=1,npar
c
c
      bma=0.5d0*dfloat(j)
      bpa=bma
c
c
c     vint : integrale de 0 a b par la methode des points de gauss.
c     ............................................................
c
      vint=0.0d0
      do 100 i=1,np2
      x=bma*xg(i)+bpa
      t=dexp(-x*x)
      t=t*hg(i)
      vint=vint+t
 100  continue
      vint=vint*bma
c
c
c     binf : integrale de b a l'infini.
c     ********************************
c
      binf=zerinf-vint
      blim=2.0d0*binf*(dfloat(j)**puist)
cCc       if(blim.lt.1.0d-14)then
cCc       if(blim.lt.1.0d-20)then
cCc       if(blim.lt.1.0d-10)then
      if(blim.lt.1.0d-20)then
                b=dfloat(j)
                go to 999
                       else
                if(j.eq.npar)then
                    write(sor,*) 'subroutine intinf : npar trop petit'
                    stop
                endif
                go to 500
      endif
c
c
 500  continue
c
c
c     l'integrale entre les bornes a et b de [u**(2z-1)*exp(-u**2) du]
c     est ensuite evaluee par la methode des points de gauss (subroutine
c     intab).
c
c
c
 999  return
      end
