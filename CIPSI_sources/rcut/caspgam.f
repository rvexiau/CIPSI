      subroutine caspgam
c
c
c     cette subroutine calcule les integrales 1/r**n dans le cas parti-
c     culier : (3), c'est-a-dire lorsque a et b sont centres en c.
c
c
c     *******************************
c     *******************************
c     evaluation de la fonction gamma :
c     *******************************
c     *******************************
c
c     n.b. : int(0,infini) signifie l'integrale entre les bornes 0 et
c            l'infini.
c
c               gamma(z)=int(0,infini)[t**(z-1)*exp(-t) dt]
c
c     ref : m. abramowitz and  i. a. stegun, "handbook of mathematical
c     ---   functions with formulas, graphs, and mathematical tables",
c           4eme edition, p.255 (1965).
c
c     en effectuant le changement de variables suivant : u**2=t, on ob-
c     tient :
c
c               gamma(z)=2*int(0,infini)[u**(2z-1)*exp(-u**2) du]
c
c
c     ..................................................................
c     on evaluera la fonction gamma entre les bornes a et l'infini  de
c     la facon suivante :
c     ..................................................................
c
c     n.b. : la borne a (rayon de coupure) est determinee de facon empi-
c            rique pour un pl donne et pour chaque atome.
c
c (a) si la puissance de u : 2z-1 est > 0 (subroutine casgam), alors on
c     aura l'egalite suivante :
c
c     2*int(a,infini)[u**(2z-1)*exp(-u**2) du]=int(0,infini)-2*int(0,a)
c
c   * int(0,infini) peut etre obtenue a l'aide de relations mathemati-
c     ques simples (cf. abramowitz p.255) :
c
c     gamma(n+1)=1*2*3...(n-1)*n=n!     avec n=entier et
c     gamma(n+1/2)=(1*3*5*7...(2n-1)/2**n)*gamma(1/2)  avec
c                                       n=nombre fractionnaire.
c
c   * int(0,a) est evaluee par la methode des points de gauss.
c
c
c
c (b) si la puissance de u : 2z-1 est < 0 ou = 0 (subroutine intinf), on
c     aura l'egalite suivante :
c
c     2*int(a,infini)[u**(2z-1)*exp(-u**2) du]=2*int(a,b)+2*int(b,inf)
c
c
c     determination de la borne b:
c     ...........................
c
c     2*int(b,inf)[exp(-u**2)/u**(k) du <
c                          (2/b**(k))*int(b,infini)[exp(-u**2)du],
c                                   avec k=(2z-1)
c     pour differentes valeurs de b, on evalue :
c     int(b,infini)=int(0,inf)[exp(-u**2)du] - int(0,b)[exp(-u**2)du],
c     la premiere integrale est egale a 1/2*(sqrt(pi)) et la seconde est
c     evaluee par la methode des points de gauss.
c     on obtient la valeur de b quand
c           (2/b**(k))*int(b,infini)[exp(-u**2)du] < 10**(-10).
c
c
c   * l'integrale entre les bornes a et b de [u**(2z-1)*exp(-u**2) du]
c     est ensuite evaluee par la methode des points de gauss (subroutine
c     intab).
c
c
c     ..................................................................
c     ..................................................................
c
c
 
      implicit  double precision(a-h,o-z)
      double precision intainf,intfgar
      integer ent,sor,pl,plder,u,casp,taba,tabb,tabpl,tabk,
     1        tabl,tabn,tabaop,tabbop,ps2a,ps2b,psab,pre1,pre2,
     1        arg,argp,plmax
      logical*1 iandj
      parameter(itypmax=3,u=20,plder=3,npar=2000000)
      integer*2 tab(15+itypmax*(2*itypmax+5))
      integer*2 tabnew
      common/tabdat/tabnew(npar)
      equivalence(tab(1),coeff)
      logical*1 iecrit
      common/ecrit/iecrit
      common/par/ijmax,npartab,plmax
      common/pmat/ pspt((itypmax+1)*(itypmax+2)*(itypmax+1)*
     1                             (itypmax+2)/4,0:4)
      common/ndebfin/ndeb(2,0:3,(itypmax+1)*(itypmax+2)/2),
     1               nfin(2,0:3,(itypmax+1)*(itypmax+2)/2)
      common/expint/itabexp((itypmax+1)*(itypmax+2)/2,
     1      (itypmax+1)*(itypmax+2)*(itypmax+1)*(itypmax+2)/4,2*itypmax)
      common/project/pl
      common/puissan/puist
      common/integr/intainf
      common/valeurz/z
      common/valpi/pi
      common/comvar/cte(2),pt
      common/puisno/no(2)
      common/norb/nf(0:itypmax)
      common/nexpres/nint,noij
      common/nocas/casp
      common/ityp/itypg,itypd
      common/vindep/ndelta,taba,tabb,tabpl,tabk,tabl,tabn,
     1              tabaop,tabbop,ndeltao
      common /viandj/iandj
      common/rayonc/rcut(0:plder)
      common/valrr/rr
      common/lecrit/ent,sor
      common/comvlab/val2a(0:u),val2b(0:u),valab(0:u)
      dimension val((itypmax+1)*(itypmax+2)*(itypmax+1)*(itypmax+2)/4)
      dimension itab(2*itypmax)
      dimension argp(10),intfgar(10)
c
c
c
c5005 format(1x,'integration numerique de la fonction gamma :',
c    1       /,1x,44('.'),/)
c
c
c
c     write(sor,5005)
c
c
c
      a=rcut(pl)*rr
c
c
      icas=0
c
c
c
      do 2500 i=1,nint
      val(i)=0.d0
 2500 continue
c
c
c
c     dans la boucle do 5000, pour chaque terme de ndeb a nfin, on
c     traite les i=1,nint integrales pour un type donne. valcte est
c     commun a toutes les nint integrales.
c
      m=ndeb(1,casp,noij)
      pre2=(m-1)*npartab
c
c     if(iecrit)then
c     write(sor,*) 'noij=',noij
c     write(sor,*) 'ndeb,nfin=',m,nfin(1,casp,noij)
c     endif
      do 5000 k=m,nfin(1,casp,noij)
c
c
      pre1=pre2+1
      pre2=pre2+npartab
      i=0
      do 3000 l=pre1,pre2
      i=i+1
      tab(i)=tabnew(l)
 3000 continue
c     if(iecrit) write(sor,*) 'tab=',(tab(ii),ii=5,npartab),coeff
c
c
      ps2a=tab(5)
      ps2b=tab(6)
      psab=tab(7)
c
      arg=tab(tabn)
      do 3200 kk=1,icas
      if(arg.eq.argp(kk)) go to 3300
 3200 continue
      icas=icas+1
      argp(icas)=arg
      vgamma=dfloat(arg-no(1))/2.0d0
      z=vgamma
      puist=(2.0d0*vgamma)-1.0d0
      if(puist.le.0.0d0)then
                call intinf
                call intab(a)
                        else
                call casgam(a)
      endif
      intfgar(icas)=intainf
      go to 3400
 3300 intainf=intfgar(kk)
c
 3400 valcte=val2a(ps2a)*val2b(ps2b)*valab(-psab)*coeff*intainf
c     if(iecrit) write(sor,*) 'intainf=',intainf
c
c
      mf1=0
      mf2=0
c
c
      do 4500 i=1,nf(itypg)
      do 4500 j=1,nf(itypd)
c
c
      mf1=mf1+1
      if(iandj .and. j.gt.i) go to 4500
c
c
c
      do 4300 n=1,2*itypmax
      itab(n)=itabexp(noij,mf1,n)
 4300 continue
c     if(iecrit) write(sor,*) 'itab=',(itab(n),n=1,2*itypmax)
c
c
      mf2=mf2+1
c
c
c     dans la boucle do 4350, on traite les deltas.
c
      do 4350 l=2,itypmax+itypd
      do 4320 m=1,l-1
      nodelta=7+(l-1)*(l-2)/2+m
      if(tab(nodelta).eq.0)then
                go to 4320
                           else
                if(itab(l).eq.itab(m))then
                                go to 4320
                                      else
                                go to 4500
                endif
      endif
 4320 continue
 4350 continue
c
c
c     la contribution du kieme terme a chacune des nint integrales est
c     mise dans le tableau val(i).
c
      val(mf2)=val(mf2)+valcte
c     if(iecrit) write(sor,*) 'mf2,val=',mf2,val(mf2)
c     if(iecrit) write(sor,*) 'valcte=',valcte
c
 4500 continue
c
c
 5000 continue
c
c
c
c     les valeurs numeriques des integrales i=1,nint sont mises dans  le
c     tableau pspt(i,j).
c
      do 8000 i=1,nint
      pspt(i,0)=cte(1)*val(i)
c     if(iecrit) write(sor,*) 'i,pspt=',i,pspt(i,0)
 8000 continue
c
c
c
      return
      end
