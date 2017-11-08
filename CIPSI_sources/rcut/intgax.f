      subroutine intgax(valaob,valaibi)
c
c
c
c     ..................................................................
c
c
c
      implicit double precision (a-h,o-z)
      double precision intj(0:10,-20:20),intj1(0:10,-20:20)
      integer ent,sor
      integer pl,plmax,xy,casp,pre1,pre2,op
      integer taba,tabb,tabaob,tabaibi,tabk,tabl,tabn,tabpl,tabaop,
     1        tabbop,tababop,ps2a,ps2b,psab,psaob,psaibi,psabop,rs,
     1        rs1,q
      logical*1 iandj
      parameter (itypmax=3,xy=20,npar=2000000)
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
      common/nocas/casp
      common/ityp/itypg,itypd
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
     1                                (itypmax+1)*(itypmax+2)/4,3)
      dimension xi(-10:10),besab(0:10),itab(2*itypmax)
      dimension nvalkln(158,2)
      dimension ndebi(0:3,2,10),nfini(0:3,2,10)
      data np1,np2/40,64/
      data ndebi/1,65,118,145,31,94,134,155,2,65,118,145,32,94,134,155,
     1           3,66,118,145,34,94,134,155,5,68,118,145,36,96,134,155,
     1           6,69,119,145,39,96,134,155,8,71,121,145,42,99,134,155,
     1           11,74,124,145,45,102,137,155,13,76,125,146,49,102,137,
     1           155,17,80,127,148,53,106,137,155,23,86,130,151,57,110,
     1           141,155/
      data nfini/1,64,117,144,31,93,133,154,2,65,117,144,33,93,133,154,
     1           4,67,117,144,35,95,133,154,5,68,118,144,38,95,133,154,
     1         7,70,120,144,41,98,133,154,10,73,123,144,44,101,136,154,
     1           12,75,124,145,48,101,136,154,16,79,126,147,52,105,136,
     1           154,22,85,129,150,56,109,140,154,30,93,133,154,64,117,
     1           144,158/
      data nvalkln /1,0,0,1,1,1,2,1,2,3,0,0,0,0,1,1,0,0,1,1,2,2,0,0,1,
     1              1,2,2,3,3,1,1,2,0,1,1,2,3,0,1,2,1,2,3,1,2,3,4,0,1,
     1              2,3,1,2,3,4,0,0,1,1,2,2,3,3,2,2,3,1,1,2,1,2,3,2,2,
     1              2,2,3,3,2,2,3,3,4,4,2,2,3,3,4,4,5,5,2,3,2,3,4,1,2,
     1              3,2,3,4,5,1,2,3,4,2,2,3,3,4,4,5,5,3,3,4,3,4,5,2,2,
     1              3,2,3,4,2,3,4,5,3,4,5,3,4,5,6,2,3,4,5,4,4,5,4,5,6,
     1              4,5,6,7,4,5,6,7,
     1             2,1,1,0,0,0,-1,0,-1,-2,-1,1,-1,1,-2,0,-1,1,-2,0,-3,
     1           -1,-1,1,-2,0,-3,-1,-4,-2,2,2,1,1,0,2,1,0,1,0,-1,0,-1,
     1            -2,2,1,0,-1,1,0,-1,-2,0,-1,-2,-3,-1,1,-2,0,-3,-1,-4,
     1            -2,1,1,0,0,0,-1,0,-1,-2,-1,1,-1,1,-2,0,-1,1,-2,0,-3,
     1             -1,-1,1,-2,0,-3,-1,-4,-2,1,0,1,0,-1,0,-1,-2,1,0,-1,
     1             -2,0,-1,-2,-3,-1,1,-2,0,-3,-1,-4,-2,0,0,-1,0,-1,-2,
     1            -1,-1,-2,-1,-2,-3,-1,-2,-3,-4,0,-1,-2,0,-1,-2,-3,-1,
     1             -2,-3,-4,-1,-1,-2,-1,-2,-3,-1,-2,-3,-4,-1,-2,-3,-4/
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
c     evaluation des integrales j(k,n) ou j(l,n) : (intj(kl,n)) par in-
c     tegration numerique avec la methode de gauss-legendre
c                                               (np2=64 points).
c     *****************************************************************
c
c
c
c     write(sor,*) 'intgax','itypg=',itypg,'itypd=',itypd
      if(casp.eq.1)then
                tababop=tabbop
                   else
                tababop=tabaop
      endif
c
c
c
c     nvalkln(158,2) : tableau contenant les valeurs des "kl" et des "n"
c                      pour chacun des pl (0,3). les ndebi(pl,casp,noij)
c                      et les nfini(pl,casp,noij) identifient le premier
c                      couple (kl,n) pour un cas particulier donne et un
c                      type noij d'integrale donne.
c
c
c     write(sor,*) 'ecriture du tableau nvalkln pour verification'
c     do 5 ll=0,3
c     write(sor,*) 'pl=ll=',ll
c     do 4 casp=1,2
c     write(sor,*) 'casp=',casp
c     do 3 ij=1,10
c     write(sor,*) 'noij=ij=',ij
c     nobori=ndebi(ll,casp,ij)
c     noborf=nfini(ll,casp,ij)
c     write(sor,*) 'nobori=',nobori,'noborf=',noborf
c     do 2 ii=nobori,noborf
c     write(sor,*)'nvalkln=',(nvalkln(ii,jj),jj=1,2)
c2    continue
c3    continue
c4    continue
c5    continue
c
c
c
      nobori=ndebi(pl,casp,noij)
      noborf=nfini(pl,casp,noij)
c
c
c     minn : valeur minimale de "n"
c     maxn : valeur maximale de "n"
c     minkl : valeur minimale de "kl"
c     maxkl : valeur maximale de "kl"
c
c
      nterm=noborf-nobori
      minn=nvalkln(nobori,2)
      maxn=minn
      if(nterm.gt.0)then
             do 10 ii=nobori+1,noborf
             if(nvalkln(ii,2).gt.maxn) maxn=nvalkln(ii,2)
             if(nvalkln(ii,2).lt.minn) minn=nvalkln(ii,2)
 10          continue
      endif
c
c
      minkl=nvalkln(nobori,1)
      maxkl=nvalkln(noborf,1)
 
c
c
c
c     integration numerique (np2=64 points) :
c     .....................................
c
c
c     x : le point courant.
c     hg(i) : le poids de gauss (weight factor).
c
c
      n1000=2
      if(dabs(bpa(1)-bpa(2)).le.1.0d-8) n1000=1
c
c
      do 30 ii=nobori,noborf
      kl=nvalkln(ii,1)
      n=nvalkln(ii,2)
      intj(kl,n)=0.0d0
      intj1(kl,n)=0.0d0
 30   continue
c
c
      do 1500 ii=1,n1000
c     write(sor,*) 'ii=',ii,'bpa(ii)=',bpa(ii),'bma(ii)=',bma(ii)
      do 1000 i=1,np2
      x=bma(ii)*xg(i)+bpa(ii)
      t=dexp((gab-x)*x-pt)
      t=t*hg(i)
      xbs=gab*x
      nbs=maxkl
      call dfbsm1
      do 100 kl=minkl,maxkl
      besab(kl)=bs(kl+1)
 100  continue
      xi(maxn)=x**(2-maxn)
      do 200 n=maxn-1,minn,-1
      xi(n)=xi(n+1)*x
 200  continue
      do 300 jj=nobori,noborf
      kl=nvalkln(jj,1)
      n=nvalkln(jj,2)
      intj(kl,n)=intj(kl,n)+besab(kl)*xi(n)*t
 300  continue
 1000 continue
c
c
      do 1250 jj=nobori,noborf
      kl=nvalkln(jj,1)
      n=nvalkln(jj,2)
      intj(kl,n)=intj(kl,n)*bma(ii)
c     write(sor,*) 'kl=',kl,'n=',n,'intj(kl,n)=',intj(kl,n)
      if(n1000.eq.1) then
                intj(kl,n)=intj(kl,n)*2.0d0
                go to 1250
      endif
      if(ii.eq.1)then
                intj1(kl,n)=intj(kl,n)
                intj(kl,n)=0.0d0
                go to 1250
      endif
      intj(kl,n)=intj(kl,n)+intj1(kl,n)
 1250 continue
c
c
 1500 continue
c
c
c
c     ------------------------------------------------------------------
c
c
c
c     ***************************************************************
c     evaluation de toutes les integrales pour un type (noij) donne :
c     ***************************************************************
c
c
c
      do 2500 i=1,nint
      do 2250 j=1,3
      val(i,j)=0.d0
 2250 continue
 2500 continue
c
c
      m=ndeb(2,casp,noij)
      pre2=(m-1)*npartab
c
      do 5000 k=m,nfin(2,casp,noij)
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
c     partie constante du terme (k) (do 5000) pour les integrales d'un
c     type noij donne :
c     ..................................................................
c
c
      psaob=tab(tabaob)
      ps2a=tab(5)
      ps2b=tab(6)
      psab=tab(7)
      rs=tab(rs1)
      q=tab(tabn)
c
      valcte=valaob(-psaob)*val2a(ps2a)*val2b(ps2b)*valab(-psab)*
     1       intj(rs,q)*coeff
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
      do 4400 tabaibi=tabaob+1,tabaob+(itypmax+itypd)
      ii=ii+1
      i1=itab(ii)
      psaibi=tab(tabaibi)
      varaibi=varaibi*valaibi(i1,psaibi)
 4400 continue
c
c
c
c     ..................................................................
c     traitement de chacun des op : x, y, z
c     ..................................................................
c
c
      do 4490 op=1,3
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
c     on evalue la contribution des aop ou des bop.
c
      psabop=tab(tababop)
      varabop=valaibi(op,psabop)
c
c
      val(mf2,op)=val(mf2,op)+valcte*varaibi*varabop
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
 7500 continue
 8000 continue
c
c
c
      return
      end
