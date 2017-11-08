      subroutine varcom
c
      implicit double precision (a-h,o-z)
      integer ent,sor,u,casp,pl,pmax2a,pmax2b,pmaxab,pmaxabi,pmaxa,
     1        pmaxb,plder,tabaob,taba,tabb,tabpl,tabk,tabl,tabn,
     1        tabaop,tabbop,rs1,dm,pli,dml,dmu
      parameter(itypmax=3,u=20,plder=3)
      include 'pshf.prm'
      common/nocas/casp
      common/ityp/itypg,itypd
      common/itypdcc/itypdc,itypgd
      common/alf/alfa,alfb
      common/rayonc/rcut(0:plder)
      common/rcrr/rc(0:plder)
      common/project/pl
      common/vindep/ndelta,taba,tabb,tabpl,tabk,tabl,tabn,
     1              tabaop,tabbop,ndeltao
      common/vintnum/pmax2a,pmax2b,pmaxab,pmaxabi,pmaxa,pmaxb
      common/dist/a(3),b(3),pa,pb
      common/distpar/ab(3),pab
      common/comvala/vala(0:u),valai(0:3,0:u)
      common/comvalb/valb(0:u),valbi(0:3,0:u)
      common/comvlab/val2a(0:u),val2b(0:u),valab(0:u)
      common/comvar/cte(2),pt
      common/vintgab/ga,gb,s
      common/puisno/no(2)
      common/borikln/bpa(0:3),bma(0:3)
      common/iklnpar/gab,tabaob,rs1
      common/iselec/ival,ivalf
      common/valpi/pi
      common/valrr/rr
      common/lecrit/ent,sor
      dimension rmin(2),rmax(2),ites(2)
      data pi4/0.1256637061435917d+2/
c
c
c
c
      valab(0)=1.d0
      val2b(0)=1.d0
      val2a(0)=1.d0
c
c
c
      rab=alfa+alfb
      rr=dsqrt(rab)
      cte(1)=dfloat(2*pl+1)*pi4/(rr**(3-no(1)))
      cte(2)=2.0d0*pi/rab
      m=casp+1
c
c
c
      go to (50,50,120,450) m
c
c
c
 50   gbb=alfb*pb
      gbb=gbb+gbb
      gb=gbb/rr
c
      valb(0)=1.d0
c
      do 55 i=0,3
      valbi(i,0)=1.0d0
 55   continue
c
      do 105 i=1,pmaxb
      valb(i)=valb(i-1)/gb
 105  continue
c
      do 115 j=1,pmaxabi
      do 110 i=1,3
      valbi(i,j)=valbi(i,j-1)*b(i)
 110  continue
 115  continue
c
c
      if(casp.eq.1) go to 400
c
 120  gaa=alfa*pa
      gaa=gaa+gaa
      ga=gaa/rr
c
      vala(0)=1.d0
c
      do 125 i=0,3
      valai(i,0)=1.0d0
 125  continue
c
      do 155 i=1,pmaxa
      vala(i)=vala(i-1)/ga
 155  continue
c
      do 165 j=1,pmaxabi
      do 160 i=1,3
      valai(i,j)=valai(i,j-1)*a(i)
 160  continue
 165  continue
c
c
c
 400  if(casp.ne.0)then
c
                go to (401,402) casp
c
c
 401            itypdc=itypd
                itypgd=itypg
                alfab=alfb
                gab=gb
                tabaob=tabb
                rs1=tabl
                go to 403
c
c
 402            itypdc=itypg
                itypgd=itypd
                alfab=alfa
                gab=ga
                tabaob=taba
                rs1=tabk
c
c
 403            pt=alfab*pab*pab
c
c
                dmu=16*itypdc
c
                do 405 i=ival,ivalf
                ites(i)=0
                if(itypgd.eq.(pl+2))then
                        dml=16*pl+16
                                    else
                        dml=16*pl
                endif
                dm=(2-no(i))*8
                dmin=dfloat(dm+dml)
                dmax=dmin+dfloat(dmu)
                if (gab*gab+dmin.lt.0.0d0)then
                        rmin(i)=0.0d0
                                          else
                        rmin(i)=dsqrt(gab*gab+dmin)
                        rmin(i)=(gab+rmin(i))*0.25d0
                endif
                if (gab*gab+dmax.lt.0.0d0)then
                        rmax(i)=0.0d0
                        ites(i)=1
                                          else
                        rmax(i)=dsqrt(gab*gab+dmax)
                        rmax(i)=(gab+rmax(i))*0.25d0
                endif
 405            continue
c
c
c
                   else
c
c
c
                pt=alfa*pa*pa+alfb*pb*pb
c
                s=ga+gb
                do 410 ii=ival,ivalf
                ites(ii)=0
                dmin=8.0d0*dfloat(2*pl+(ii+1)-no(ii))
                dmax=dmin+dfloat(16*(itypg+itypd))
                if (s*s+dmin.lt.0.0d0)then
                        rmin(ii)=0.0d0
                                      else
                        rmin(ii)=dsqrt(s*s+dmin)
                        rmin(ii)=(rmin(ii)+s)*0.25d0
                endif
                if (s*s+dmax.lt.0.0d0)then
                        rmax(ii)=0.0d0
                        ites(ii)=1
                                      else
                        rmax(ii)=dsqrt(s*s+dmax)
                        rmax(ii)=(rmax(ii)+s)*0.25d0
                endif
 410            continue
c
c
c
      endif
c
      do 415 i=ival,ivalf
      pli=pl+1
      if(pli.gt.plder) pli=plder
      do 412 iii=pl,pli
      rc(iii)=rcut(iii)*rr
 412  continue
c     write(sor,*) 'i=',i,'rmin(i)=',rmin(i)
      if(rmin(i).gt.6.0d0) rc(pl)=dmax1(rc(pl),rmin(i)-6.0d0)
c     write(sor,*) 'pl=',pl,'rc(pl)=',rc(pl)
      if (ites(i).eq.0)then
                bpa(i-1)=(rc(pl)+rmax(i)+6.0d0)*0.5d0
                bma(i-1)=(rmax(i)+6.0d0-rc(pl))*0.5d0
                       else
                bpa(i-1)=(2.0d0*rc(pl)+6.0d0)*0.5d0
                bma(i-1)=3.0d0
      endif
      if(ivalf.eq.1) go to 415
      if(i.eq.2)then
             if(pl.eq.3)then
                        bpa(2)=bpa(1)
                        bma(2)=bma(1)
                else
		if(rcut(pl+1).eq.rcut(pl) )then
                        bpa(2)=bpa(1)
                        bma(2)=bma(1)
                                      else
             if(rmin(i).gt.6.0d0) rc(pl+1)=dmax1(rc(pl+1),rmin(i)-6.0d0)
                        bpa(2)=(rc(pl+1)+rmax(i)+6.0d0)*0.5d0
                        bma(2)=(rmax(i)+6.0d0-rc(pl+1))*0.5d0
                end if
                endif
      endif
 415  continue
c
c
c
 450  do 500 j=1,pmaxab
      valab(j)=valab(j-1)/rab
 500  continue
c
      do 550 j=1,pmax2b
      val2b(j)=val2b(j-1)*2*alfb
 550  continue
c
      do 600 j=1,pmax2a
      val2a(j)=val2a(j-1)*2*alfa
 600  continue
c
c
c
      return
      end
