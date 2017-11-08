      subroutine caspgax
c
c
c     ..................................................................
c
c
c
 
      implicit  double precision(a-h,o-z)
      double precision intainf,int12,intfgar
      integer ent,sor,pl,plder,u,casp,taba,tabb,tabpl,tabk,
     1        tabl,tabn,tabaop,tabbop,ps2a,ps2b,psab,pre1,pre2,
     1        arg,argp,op,plmax
      logical*1 iandj
      parameter(itypmax=3,plder=3,u=20,npar=2000000)
      integer*2 tab(15+itypmax*(2*itypmax+5))
      integer*2 tabnew
      common/tabdat/tabnew(npar)
      equivalence(tab(1),coeff)
      common/par/ijmax,npartab,plmax
      common/pmat/ pspt((itypmax+1)*(itypmax+2)*(itypmax+1)*
     1                             (itypmax+2)/4,0:4)
      common/ndebfin/ndeb(2,0:3,(itypmax+1)*(itypmax+2)/2),
     1               nfin(2,0:3,(itypmax+1)*(itypmax+2)/2)
      common/expint/itabexp((itypmax+1)*(itypmax+2)/2,
     1      (itypmax+1)*(itypmax+2)*(itypmax+1)*(itypmax+2)/4,2*itypmax)
      common/puissan/puist
      common/integr/intainf
      common/valeurz/z
      common/comvar/cte(2),pt
      common/norb/nf(0:itypmax)
      common/project/pl
      common/nexpres/nint,noij
      common/nocas/casp
      common/ityp/itypg,itypd
      common/rayonc/rcut(0:plder)
      common/valrr/rr
      common/vindep/ndelta,taba,tabb,tabpl,tabk,tabl,tabn,
     1              tabaop,tabbop,ndeltao
      common /viandj/iandj
      common/lecrit/ent,sor
      common/comvlab/val2a(0:u),val2b(0:u),valab(0:u)
      dimension val((itypmax+1)*(itypmax+2)*(itypmax+1)*(itypmax+2)/4,3)
      dimension itab(2*itypmax)
      dimension argp(10),int12(2),intfgar(10),a(2)
c
c
c
c5005 format(1x,'integration numerique de la fonction gamma :',
c    1       /,1x,44('.'),/)
c
c
c
      do 100 i=1,2
      int12(i)=0.0d0
 100  continue
c
c
      a(1)=rcut(pl)*rr
      if(pl.eq.3)then
                a(2)=a(1)
                 else
                a(2)=rcut(pl+1)*rr
      endif
c
c
      n1000=2
      if(dabs(a(1)-a(2)).lt.1.0d-8)n1000=1
c
c
c
      icas=0
c
c
c
      do 2500 i=1,nint
      do 2000 j=1,3
      val(i,j)=0.d0
 2000 continue
 2500 continue
c
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
      ps2a=tab(5)
      ps2b=tab(6)
      psab=tab(7)
c
c
      arg=tab(tabn)
      do 3200 kk=1,icas
      if(arg.eq.argp(kk)) go to 3300
 3200 continue
      icas=icas+1
      argp(icas)=arg
      vgamma=dfloat(arg)/2.0d0
      z=vgamma
      puist=(2.0d0*vgamma)-1.0d0
      if(puist.le.0.0d0)then
                call intinf
                do 3250 ii=1,n1000
                call intab(a(ii))
                int12(ii)=intainf
 3250           continue
                        else
                do 3275 ii=1,n1000
                call casgam(a(ii))
                int12(ii)=intainf
 3275           continue
      endif
      if(n1000.eq.1)then
                intfgar(icas)=int12(1)*2.0d0
                    else
                intfgar(icas)=int12(1)+int12(2)
      endif
      intainf=intfgar(icas)
      go to 3400
 3300 intainf=intfgar(kk)
c
c
c
 3400 valcte=val2a(ps2a)*val2b(ps2b)*valab(-psab)*coeff*intainf
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
      do 4490 op=1,3
c
c     on traite les deltas op. des que l'on trouve un delta qui annule
c     la contribution du terme pour un op specifique, on passe a un
c     autre op.
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
      val(mf2,op)=val(mf2,op)+valcte
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
c     les valeurs numeriques des integrales i=1,nint sont mises dans  le
c     tableau pspt(i,j).
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
