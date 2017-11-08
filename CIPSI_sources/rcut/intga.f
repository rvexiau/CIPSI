      subroutine intga(valaob,valaibi)
c
c
c     cette subroutine calcule les integrales 1/r**n dans les cas
c     particuliers : (1) a est centre en c et (2) b est centre en c.
c     dans la premiere partie, on evalue tout d'abord les j(k,n) ou
c     j(l,n) par integration numerique avec la methode de gauss-
c     legendre. ensuite, dans la seconde partie, on calcule toutes
c     les integrales (i=1,nint) pour un type donne.
c
c     la valeur numerique de chaque integrale (i), pour un type donne,
c     est mise dans pspt(i,0).
c
c
c *** pour des definitions et explications detaillees, cf. subroutine
c     intgab.
c
c     ..................................................................
c
c
c
      implicit double precision (a-h,o-z)
      double precision intj(0:10,-20:20)
      integer ent,sor
      integer pl,xy,casp,plder,pre1,pre2
      integer taba,tabb,tabaob,tabaibi,tabk,tabl,tabn,tabpl,tabaop,
     1        tabbop,ps2a,ps2b,psab,psaob,psaibi,rs,rs1,q
      logical*1 iandj
      parameter (itypmax=3,xy=20,plder=3,npar=2000000)
      integer*2 tab(15+itypmax*(2*itypmax+5))
      integer*2 tabnew
      common/tabdat/tabnew(npar)
      equivalence (tab(1),coeff)
      common/uncp/xbs,bs(10),nbs
      common/numeriq/xg(64),hg(64)
      common/pmat/ pspt((itypmax+1)*(itypmax+2)*(itypmax+1)*
     1                              (itypmax+2)/4,0:4)
      common/ndebfin/ndeb(2,0:3,(itypmax+1)*(itypmax+2)/2),
     1               nfin(2,0:3,(itypmax+1)*(itypmax+2)/2)
      common/expint/itabexp((itypmax+1)*(itypmax+2)/2,
     1      (itypmax+1)*(itypmax+2)*(itypmax+1)*(itypmax+2)/4,2*itypmax)
      common/comvlab/val2a(0:xy),val2b(0:xy),valab(0:xy)
      common/comvar/cte(2),pt
      common/par/ijmax,npartab,plmax
      common/puisno/no(2)
      common/nocas/casp
      common/ityp/itypg,itypd
      common/itypdcc/itypdc,itypgd
      common/norb/nf(0:itypmax)
      common/nexpres/nint,noij
      common/vindep/ndelta,taba,tabb,tabpl,tabk,tabl,tabn,
     1              tabaop,tabbop,ndeltao
      common/project/pl
      common /viandj/iandj
      common/borikln/bpa(0:3),bma(0:3)
      common/iklnpar/gab,tabaob,rs1
      common/lecrit/ent,sor
      dimension valaob(0:xy),valaibi(0:3,0:xy)
      dimension val((itypmax+1)*(itypmax+2)*
     1                                (itypmax+1)*(itypmax+2)/4)
      dimension xi(-10:10),itab(2*itypmax)
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
c     evaluation des integrales j(k,n) ou j(l,n) : (intj(kl,n)) par in-
c     tegration numerique avec la methode de gauss-legendre
c                                               (np2=64 points).
c     *****************************************************************
c
c
c      definitions de variables generales utilisees soit dans le cas
c      particulier (1), soit dans le cas particulier (2) :
c      ................................................................
c
c
c      klmax : valeur maximale de k ou de l dans les integrales j(k,n)
c              et j(l,n).
c
c      nmax : valeur maximale de n dans les integrales j(k,n) ou j(l,n).
c
c
c
c
      no1=no(1)
c
c
c
      klmax=pl+itypdc
      if(itypgd.eq.(pl+2))then
                nmax=no1-klmax-2
                n1=no1-2
                          else
                nmax=no1-klmax
                n1=no1
      endif
c
c
c
c     definitions des bornes d'integration pour les j(k,n) ou j(l,n) :
c     ..............................................................
c
c
c
c     integration numerique (np2=64 points) :
c     .....................................
c
c
c     notons que l'on calcule dans la boucle do 1000 toutes les inte-
c     grales j(k,n) et j(l,n) necessaires pour l'evaluation numeri-
c     que d'un type d'integrale donne.
c
c     x : le point courant.
c     hg(i) : le poids de gauss (weight factor).
c
c
c
      do 30 kl=pl,klmax
      n=n1-kl
      intj(kl,n)=0.0d0
 30   continue
c
c
      do 1000 i=1,np2
      x=bma(0)*xg(i)+bpa(0)
      t=dexp((gab-x)*x-pt)
      t=t*hg(i)
      xbs=gab*x
      nbs=pl+itypdc
      call dfbsm1
      xi(no1)=x**(2-no1)
      do 200 n=no1-1,nmax,-1
      xi(n)=xi(n+1)*x
 200  continue
      n=n1-pl+1
      do 300 kl=pl,klmax
      n=n-1
      intj(kl,n)=intj(kl,n)+bs(kl+1)*xi(n)*t
 300  continue
 1000 continue
c
c
      do 1300 kl=pl,klmax
      n=n1-kl
      intj(kl,n)=intj(kl,n)*bma(0)
 1300 continue
c
c
c
c     ------------------------------------------------------------------
c
c
c
c     ********************************************************
c     evaluation de toutes les integrales pour un type donne :
c     ********************************************************
c
c
c     valaob(0:xy) : valeur des parametres a ou b du tableau tab(npartab).
c     valaibi(0:3,0:xy) : valeurs des ai et des bi selon le cas pour des
c                         puissances de 0 a une valeur maximale (2eme in-
c                         dice). chacune de ces puissances est couplee a :
c                         0 pour s, 1 pour x, 2 pour y et 3 pour z, ces
c                         derniers etant definis par le premier indice de
c                         ce tableau.
c
c
c
      do 2500 i=1,nint
      val(i)=0.d0
 2500 continue
c
c
      m=ndeb(1,casp,noij)
      pre2=(m-1)*npartab
c
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
c
c
c     ..................................................................
c     partie constante du terme (k) pour chacun des termes d'un type
c     d'integrale donne :
c     ..................................................................
c
c
      psaob=tab(tabaob)
      ps2a=tab(5)
      ps2b=tab(6)
      psab=tab(7)
      rs=tab(rs1)
      q=tab(tabn)+no1
c
      valcte=valaob(-psaob)*val2a(ps2a)*val2b(ps2b)*valab(-psab)*
     1       intj(rs,q)*coeff
c
c
c     ..................................................................
c     evaluation de la contribution des ai et des bi pour chacune des
c     integrales (i=1,nint) d'un type donne :
c     ..................................................................
c
c
c     cf. subroutine intgab
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
      do 4300 n=1,2*itypmax
      itab(n)=itabexp(noij,mf1,n)
 4300 continue
c
c
      mf2=mf2+1
c
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
      ii=0
      varaibi=1.d0
      do 4400 tabaibi=tabaob+1,tabaob+2*itypmax
      ii=ii+1
      i1=itab(ii)
      psaibi=tab(tabaibi)
      varaibi=varaibi*valaibi(i1,psaibi)
 4400 continue
c
      val(mf2)=val(mf2)+valcte*varaibi
c
 4500 continue
c
c
 5000 continue
c
c
c
c     pspt(i,j) pour i=1,nint : valeur numerique de chaque integrale
c                             (i) pour un type donne.
c
c
c     write(sor,5003)
      do 8000 i=1,nint
      pspt(i,0)=cte(1)*val(i)
 8000 continue
c
c
c
      return
      end
