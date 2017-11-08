      subroutine intgab
c
c
c     cette subroutine calcule les integrales 1/r**n dans le cas ge-
c     neral. dans la premiere partie, on evalue tout d'abord les
c     i(k,l,n) par integration numerique avec la methode de gauss-
c     legendre. ensuite, dans la seconde partie, on calcule toutes
c     les integrales (i=1,nint) pour un type donne.
c
c     la valeur numerique de chaque integrale (i), pour un type donne
c     est mise dans pspt(i).
c
c
c     .................................................................
c
c
      implicit  double precision(a-h,o-z)
      integer ent,sor,u,casp,pl,r,ss,q,taba,tabb,tabpl,tabk,tabl,tabn,
     1        psa,psb,ps2a,ps2b,psab,psai,psbi,tabai,tabbi,tabbop,
     1        tabaop,plder,pre1,pre2,plmax
      double precision int
      logical*1 iandj
      parameter(itypmax=3,u=20,plder=3,npar=2000000)
      common/uncp/xbs,bs(10),nbs
      common/numeriq/xg(64),hg(64)
      common/expint/itabexp((itypmax+1)*(itypmax+2)/2,
     1      (itypmax+1)*(itypmax+2)*(itypmax+1)*(itypmax+2)/4,2*itypmax)
      common/pmat/ pspt((itypmax+1)*(itypmax+2)*(itypmax+1)*
     1                             (itypmax+2)/4,0:4)
      common/nocas/casp
      common/ityp/itypg,itypd
      common/norb/nf(0:itypmax)
      common/nexpres/nint,noij
      common/vindep/ndelta,taba,tabb,tabpl,tabk,tabl,tabn,
     1              tabaop,tabbop,ndeltao
      common/ndebfin/ndeb(2,0:3,(itypmax+1)*(itypmax+2)/2),
     1               nfin(2,0:3,(itypmax+1)*(itypmax+2)/2)
      common/comvar/cte(2),pt
      common/par/ijmax,npartab,plmax
      common/project/pl
      common/puisno/no(2)
      common /viandj/iandj
      common/comvala/vala(0:u),valai(0:3,0:u)
      common/comvalb/valb(0:u),valbi(0:3,0:u)
      common/comvlab/val2a(0:u),val2b(0:u),valab(0:u)
      common/vintgab/ga,gb,s
      common/borikln/bpa(0:3),bma(0:3)
      common/lecrit/ent,sor
      dimension int(0:10,0:10,-20:20)
      dimension xi(-10:10),besa(0:10),besb(0:10)
      dimension val((itypmax+1)*(itypmax+2)*(itypmax+1)*(itypmax+2)/4)
      dimension itab(2*itypmax)
      integer*2 tabnew
      common/tabdat/tabnew(npar)
      common/varvpl/vpl(0:3)
      integer*2 tab(15+itypmax*(2*itypmax+5))
      equivalence (tab(1),coeff)
      data np1,np2/40,64/
c
c
c
c5001 format(1x,///,1x,'resultats :',/,1x,11('*'),//)
c     write(sor,5001)
 5003 format(1x,///,1x,'pspt(i,j) :',/,1x,9('*'),//)
c
c
c
c     *****************************************************************
c     evaluation des integrales i(k,l,n) : (int(k,l,n)) par integration
c     numerique avec la methode de gauss-legendre (np2=64 points).
c     *****************************************************************
c
c
c
c
      no1=no(1)
c
c
c
c     integration numerique (np2=64 points) :
c     .....................................
c
c
c     notons que l'on calcule dans la boucle do 1000 toutes les inte-
c     grales i(k,l,n) necessaires pour l'evaluation numerique d'un
c     type d'integrale donne.
c
c     x : le point courant
c     hg(i) : le poids de gauss (weight factor)
c
c
c
      do 50 k=pl,pl+itypg
      do 40 l=pl,pl+itypd
      n=2*pl+no1-k-l
      int(k,l,n)=0.0d0
40    continue
50    continue
c
c
      do 1000 i=1,np2
      x=bma(0)*xg(i)+bpa(0)
      t=dexp((s-x)*x-pt)
      t=t*hg(i)
      xbs=ga*x
      nbs=pl+itypg
      call dfbsm1
      do 100 k=pl,pl+itypg
      besa(k)=bs(k+1)
 100  continue
      xbs=gb*x
      nbs=pl+itypd
      call dfbsm1
      do 200 l=pl,pl+itypd
      besb(l)=bs(l+1)
 200  continue
      xi(no1)=x**(2-no1)
      do 300 n=no1-1,no1-(itypg+itypd),-1
      xi(n)=xi(n+1)*x
 300  continue
      do 500 k=pl,pl+itypg
      do 400 l=pl,pl+itypd
      n=2*pl+no1-k-l
      int(k,l,n)=int(k,l,n)+besa(k)*besb(l)*xi(n)*t
 400  continue
 500  continue
 1000 continue
c
c
      do 1500 k=pl,pl+itypg
      do 1400 l=pl,pl+itypd
      n=2*pl+no1-k-l
      int(k,l,n)=int(k,l,n)*bma(0)
 1400 continue
 1500 continue
c
c
c
c......................................................................
c
c
c
c     ********************************************************
c     evaluation de toutes les integrales pour un type donne :
c     ********************************************************
c
c
c
      do 2500 i=1,nint
      val(i)=0.d0
 2500 continue
c
c
c
c     pour un type d'integrale donne correspondant au numero noij, on
c     repere les ndeb(noij,casp) et nfin(noij,casp). on lit ensuite
c     chacun des termes et l'on evalue la partie constante du terme
c     (valcte). pour chacune des integrales d'un type donne, dans
c     un deuxieme temps on calcule la partie qui depend des ai et
c     des bi (varai et varbi). le produit valcte*varai*varbi donne
c     ainsi la contribution du terme a chacune des integrales (1,nint)
c     du type en question. on fait une sommation sur tous les termes
c     pour chacune des integrales du type traite.
c
c
      m=ndeb(1,casp,noij)
      pre2=(m-1)*npartab
c
c
      do 5000 k=m,nfin(1,casp,noij)
c
      pre1=pre2+1
      pre2=pre2+npartab
      i=0
      do 3000 l=pre1,pre2
      i=i+1
      tab(i)=tabnew(l)
 3000 continue
c
c
c     partie constante du terme pour chacun des termes d'un type d'inte-
c     grale donne :
c     ..................................................................
c
      psa=tab(taba)
      psb=tab(tabb)
      ps2a=tab(5)
      ps2b=tab(6)
      psab=tab(7)
      nopl=tab(tabpl)
      r=tab(tabk)
      ss=tab(tabl)
      q=tab(tabn)+no1
c
      valcte=vala(-psa)*valb(-psb)*val2a(ps2a)*val2b(ps2b)*
     1       valab(-psab)*vpl(nopl)*int(r,ss,q)*coeff
c
c
c     ..................................................................
c     evaluation de la contribution des ai et des bi pour chacune des
c     integrales (i=1,nint) d'un type donne :
c     ..................................................................
c
c
c
      mf1=0
      mf2=0
c
c
      do 4500 i=1,nf(itypg)
      do 4500 j=1,nf(itypd)
c
      mf1=mf1+1
      if(iandj .and. j.gt.i) go to 4500
c
      do 4300 n=1,2*itypmax
      itab(n)=itabexp(noij,mf1,n)
 4300 continue
c
      mf2=mf2+1
c
c *** analyse des deltas. des que l'on trouve un delta qui annulle le
c     terme pour une integrale (i), on passe a la suivante.
c
c     nodelta : position du delta dans le tableau tab(npartab).
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
c     l'expression numerique qui decrit chacune des integrales (1,nint)
c     est mise dans le tableau itab(ii) ou itab(jj). ii et jj sont des
c     compteurs qui couplent la position de chaque nombre du tableau
c     itab avec les ai ou bi selon le cas. on peut ensuite aller cher-
c     cher dans le tableau valai(i1,psai), les valeurs des ai ou dans
c     valbi(i1,psbi), les valeurs des bi.
c
c
      ii=0
      varai=1.d0
      do 4400 tabai=taba+1,taba+2*itypmax
      ii=ii+1
      i1=itab(ii)
      psai=tab(tabai)
      varai=varai*valai(i1,psai)
 4400 continue
c
      jj=0
      varbi=1.d0
      do 4450 tabbi=tabb+1,tabb+2*itypmax
      jj=jj+1
      j1=itab(jj)
      psbi=tab(tabbi)
      varbi=varbi*valbi(j1,psbi)
 4450 continue
c
c
      val(mf2)=val(mf2)+valcte*varai*varbi
c
 4500 continue
c
c
 5000 continue
c
c
c
c     pspt(i,j) pour i=1,nint : valeur numerique de chaque integrale
c                               (i) pour un type donne. l'indice j de-
c                               finit si l'on traite : 1/r**4 ou
c                                                x,y,z/r**3.
c
      do 8000 i=1,nint
      pspt(i,0)=cte(1)*val(i)
 8000 continue
c
c
c
      return
      end
