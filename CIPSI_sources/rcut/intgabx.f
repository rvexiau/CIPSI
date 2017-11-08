      subroutine intgabx
c
c
c     .................................................................
c
c
      implicit  double precision(a-h,o-z)
      integer ent,sor,u,casp,pl,r,ss,q,taba,tabb,tabpl,tabk,tabl,tabn,
     1        psa,psb,ps2a,ps2b,psab,psai,psbi,psaop,psbop,tabai,tabbi,
     1        tabbop,tabaop,pre1,pre2,plmax,op
      double precision  int(0:10,0:10,-20:20),int1(0:10,0:10,-20:20)
      double precision  intest(0:10,0:10,-20:20)
      logical*1 iandj,iecrit
      parameter(itypmax=3,u=20,npar=2000000)
      common/ecrit/iecrit
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
      dimension xi(-10:10),besa(0:10),besb(0:10)
      dimension val((itypmax+1)*(itypmax+2)*(itypmax+1)*(itypmax+2)/4,3)
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
 5001 format(1x,///,1x,'resultats :',/,1x,11('*'),//)
c     write(sor,5001)
 5003 format(1x,///,1x,'pspt(i,j) :',/,1x,9('*'),//)
c
c
c
c     *****************************************************************
c     *****************************************************************
c     evaluation des integrales i(k,l,n) : (int(k,l,n)) par integration
c     numerique avec la methode de gauss-legendre (np2=64 points).
c     *****************************************************************
c     *****************************************************************
c
c
c
c     integration numerique (np2=64 points) :
c     .....................................
c
c
c     notons que l'on calcule dans la boucle do 1500 les integrales
c     i(k,l,n) communes a deux ("l") avec un rayon de coupure diffe-
c     rent pour chacun de ces "l".
c
c     ex : l=0      i(1,0,2), i(0,1,2)
c                      avec bpa(1) et bma(1) pour l=0
c
c          l=1      i(1,0,2), i(0,1,2)
c                      avec bpa(2) et bma(2) pour l=1
c
c     si bpa et bma sont egaux pour les deux "l" distincts parce que
c     les rayons de coupure assignes a ces "l" sont par hazard egaux,
c     on mulitiple chacune des int(k,l,n) par 2. il en est de meme dans
c     le cas ou "l"=3 afin de tronquer la serie correctement.
c
c
c
c     x : le point courant
c     hg(i) : le poids de gauss (weight factor)
c
c
c     write(sor,*) 'intgabx','itypg=',itypg,'itypd=',itypd
      n1000=2
      if(dabs(bpa(1)-bpa(2)).le.1.0d-8) n1000=1
c
c
c
      no2=no(2)
c
c
c
      do 25 k=pl,pl+itypg+1
      do 20 l=pl,pl+itypd+1
      n=2*pl+no2+1-k-l
      int(k,l,n)=0.0d0
      int1(k,l,n)=0.0d0
 20   continue
 25   continue
c
c
c
c     write(sor,*) 'n1000=',n1000
      do 1500 ii=1,n1000
c     write(sor,*) 'ii=',ii,'bpa(ii)=',bpa(ii),'bma(ii)=',bma(ii)
      do 1000 i=1,np2
      x=bma(ii)*xg(i)+bpa(ii)
      t=dexp((s-x)*x-pt)
      t=t*hg(i)
      xbs=ga*x
      nbs=pl+itypg+1
      call dfbsm1
      do 100 k=pl,pl+itypg+1
      besa(k)=bs(k+1)
 100  continue
      xbs=gb*x
      nbs=pl+itypd+1
      call dfbsm1
      do 200 l=pl,pl+itypd+1
      besb(l)=bs(l+1)
 200  continue
      xi(no2+1)=x**(1-no2)
      do 300 n=no2,no2-1-(itypg+itypd),-1
      xi(n)=xi(n+1)*x
 300  continue
      do 500 k=pl,pl+itypg+1
      do 400 l=pl,pl+itypd+1
      n=2*pl+no2+1-k-l
      int(k,l,n)=int(k,l,n)+besa(k)*besb(l)*xi(n)*t
      if(i.lt.3 .or. i.gt.62)then
                        intest(k,l,n)=besa(k)*besb(l)*xi(n)*t
c     write(sor,*) 'i=',i,'k=',k,'l=',l,'n=',n,'int=',intest(k,l,n)
      endif
 400  continue
 500  continue
 1000 continue
c
c
      do 1250 k=pl,pl+itypg+1
      do 1200 l=pl,pl+itypd+1
      n=2*pl+no2+1-k-l
      int(k,l,n)=int(k,l,n)*bma(ii)
c     write(sor,*) 'k=',k,'l=',l,'n=',n,'int(k,l,n)=',int(k,l,n)
      if(n1000.eq.1) then
                int(k,l,n)=int(k,l,n)*2.0d0
                go to 1200
      endif
      if(ii.eq.1)then
                int1(k,l,n)=int(k,l,n)
                int(k,l,n)=0.0d0
                go to 1200
      endif
      int(k,l,n)=int(k,l,n)+int1(k,l,n)
 1200 continue
 1250 continue
c
c
 1500 continue
c
c
c
c......................................................................
c
c
c
c     ***************************************************************
c     ***************************************************************
c     evaluation de toutes les integrales pour un type (noij) donne :
c     ***************************************************************
c     ***************************************************************
c
c
c
      do 2500 i=1,nint
      do 2450 j=1,3
      val(i,j)=0.0d0
 2450 continue
 2500 continue
c
c
c
      m=ndeb(2,casp,noij)
      pre2=(m-1)*npartab
c
c
c     if(iecrit) write(sor,*) 'ndeb,nfin=',m,nfin(2,casp,noij)
      do 5000 k=m,nfin(2,casp,noij)
c     if(iecrit) write(sor,*)'k=',k
c
      pre1=pre2+1
      pre2=pre2+npartab
      i=0
      do 3000 l=pre1,pre2
      i=i+1
      tab(i)=tabnew(l)
 3000 continue
c
c     if(iecrit) write(sor,*) (tab(hh),hh=5,npartab),coeff
c
c
c     ..................................................................
c     partie constante du terme (k) (do 5000) pour les integrales d'un
c     type noij donne :
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
      q=tab(tabn)
c
      valcte=vala(-psa)*valb(-psb)*val2a(ps2a)*val2b(ps2b)*
     1       valab(-psab)*vpl(nopl)*int(r,ss,q)*coeff
c
c
c     ..................................................................
c     evaluation de la contribution des ai et des bi pour chacune des
c     integrales formees dans la boucle do 4500 par les couples (i,j).
c     notons que le test : if(iandj .and. j.gt.i) evite de calculer deux
c     fois les memes integrales lorsque itypg et itypd sont egaux.
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
c     if(iecrit) write(sor,*) 'itab=',(itab(n),n=1,2*itypmax)
c
      mf2=mf2+1
c
c *** analyse des deltas. des que l'on trouve un delta qui annulle le
c     terme pour une integrale formee par do 4500, on passe a la sui-
c     te.
c
c
c     ex :         si itypmax=3
c
c     dans le tableau tab(npartab), les positions des deltas sont les
c     suivantes :
c
c      8   9  10  11  12  13  14  15  16  17  18  19  20  21  22
c     ji  ki  kj  li  lj  lk  mi  mj  mk  ml  ni  nj  nk  nl  nm
c
c     ces positions peuvent etre calculees comme l'adresse des cases de
c     la demie-matrice suivante (sans la diagonale) :
c
c                  1  2  3  4  5  6
c                  i  j  k  l  m  n
c          1  i   ii
c          2  j   ji jj
c          3  k   ki kj kk
c          4  l   li lj lk ll
c          5  m   mi mj mk ml mm
c          6  n   ni nj nk nl nm nn
c
c
c     avec nodelta=7+(l-1)*(l-2)/2+m     ou l=lignes et m=colonnes.
c
c     le tableau itab contenant les 2*itypmax nombres qui decrivent
c     l'integrale precise a calculer a l'exemple de 110 100 pour une
c     ij/l (d(xx)/p(x)), il suffit de comparer itab(l) et itab(m).
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
c
c
      ii=0
      varai=1.d0
      do 4400 tabai=taba+1,taba+(itypmax+itypd)
      ii=ii+1
      i1=itab(ii)
      psai=tab(tabai)
      varai=varai*valai(i1,psai)
 4400 continue
c
      jj=0
      varbi=1.d0
      do 4450 tabbi=tabb+1,tabb+(itypmax+itypd)
      jj=jj+1
      j1=itab(jj)
      psbi=tab(tabbi)
      varbi=varbi*valbi(j1,psbi)
 4450 continue
c
c
c     ..................................................................
c     traitement de chacun des op : x, y, z
c     ..................................................................
c
c
      do 4490 op=1,3
c
c     on traite tout d'abord les deltas op. des que l'on trouve un delta
c     qui annule la contribution du terme pour un op specifique, on pas-
c     se a un autre op.
c
      do 4455 ll=ndeltao+1,ndeltao+(itypmax+itypd)
      if(tab(ll).eq.0)then
                go to 4455
                      else
                if(op.eq.itab(ll-ndeltao))then
                                go to 4455
                                          else
                                go to 4490
                endif
      endif
 4455 continue
c
c
c     on evalue la contribution des aop et des bop.
c
      psaop=tab(tabaop)
      psbop=tab(tabbop)
      varaop=valai(op,psaop)
      varbop=valbi(op,psbop)
c
c
      val(mf2,op)=val(mf2,op)+valcte*varai*varbi*varaop*varbop
c     if(iecrit) write(sor,*)'mf2,op=',mf2,op,'term=',valcte*varai*
c    1                         varbi*varaop*varbop
c
c
 4490 continue
c
c
 4500 continue
c
c
 5000 continue
c
c
c
      do 8000 i=1,nint
      do 7500 j=1,3
      pspt(i,j)=cte(2)*val(i,j)
c     if(iecrit) write(sor,*) 'i=',i,'j=',j,'pspt=',pspt(i,j)
 7500 continue
 8000 continue
c
c
c
      return
      end
