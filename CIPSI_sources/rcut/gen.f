      subroutine gen
c
c
c     ------------------------------------------------------------------
c
c
c
      implicit integer (a-z)
      double precision xg,hg,pi
      parameter (itypmax=3)
      common/expint/itabexp((itypmax+1)*(itypmax+2)/2,
     1      (itypmax+1)*(itypmax+2)*(itypmax+1)*(itypmax+2)/4,2*itypmax)
      common/norb/nf(0:itypmax)
      common/vindep/ndelta,taba,tabb,tabpl,tabk,tabl,tabn,
     1              tabaop,tabbop,ndeltao
      common/numeriq/xg(64),hg(64)
      common/valpi/pi
      dimension e(0:itypmax),per(35,itypmax),itypab(2)
c
c     ..................................................................
c     ***********                                                      .
c     * warning *                                                      .
c     ***********                                                      .
c                                                                      .
c     lorsque itypmax est plus grand que 4 (orbitale g), il faut aug-  .
c     menter le premier indice du tableau per.                         .
c     ..................................................................
c
      dimension niexp(2*itypmax)
      dimension tabor (20)
      data tabor /1,2,3,4,5,7,10,6,8,9,11,14,20,12,15,13,17,18,19,16/
c
c
c
      call psepin
c
c
c
c     nf(i) : tableau contenant le nombre d'orbitales pour un type
c             donne.
c             i=0 (s), i=1 (p), i=2 (d), i=3 (f), ...
c             par exemple dans le cas d'une f, nf(3)=10.
c
c
      nf(0)=1
      do 2 i=1,itypmax
      nf(i)=(i+1)*(i+2)/2
 2    continue
c
      do i=1,35
      do j=1,itypmax
      per(i,j)=0
      end do
      end do
c
c
      ndelta=itypmax*(2*itypmax+1)
      taba=12+ndelta
      tabb=14+ndelta+2*itypmax
      tabpl=11+ndelta
      tabk=8+ndelta
      tabl=9+ndelta
      tabn=10+ndelta
c
c
c
c     tabaop et tabbop sont les numeros des cases du tableau tab qui
c     correspondent a aop et bop respectivement.
c     ndeltao+1 : numero du premier delta op dans le tableau tab.
c         ndeltao+1=7+((2*itypmax*2*itypmax)-(2*itypmax))/2+1
c
c
      tabaop=taba+2*itypmax+1
      tabbop=tabb+2*itypmax+1
      ndeltao=7+((4*itypmax*itypmax)-(2*itypmax))/2
c
c
c
      pi=4.0d0*datan(1.0d0)
c
c
c
c     on genere les expressions numeriques decrivant les differents ty-
c     pes d'orbitales de (0 a itypmax). les expressions des integrales
c     <ita/pl/itb> sont ensuite construites.
c     l'expression d'une ieme integrale  est formee de 2*itypmax nombres.
c     par convention, on pose s=0, x=1, y=2 et z=3.
c     on a par exemple pour l'integrale : <d/s> avec (d=xx) et itypmax=3,
c     l'expression suivante : 1 1 0 0 0 0 .
c
c
c     le resultat est mis dans itabexp.
c
c
c     .................................................................
c *** dans un premier temps, on genere les expressions numeriques decri-
c     vant le type d'orbitale (de l'orbitale s jusqu'a celle qui corres-
c     pond au itypmax defini par l'enonce parameter).
c     ..................................................................
c
c
      do 5 j=1,itypmax
      e(j)=0
      per(1,j)=0
 5    continue
c
c
c
      cont=1
      e(0)=3
      ityp=0
c
c
c
 10   ityp=ityp+1
      if(ityp.gt.itypmax) go to 100
c
      do 15 j=1,ityp
      e(j)=0
 15   continue
c
c
c     i : numero de la boucle.
c     e(i) : valeur de l'indice de la boucle i.
c     per(cont,j) : tableau contenant les valeurs numeriques qui de-
c                   crivent l'expression d'une orbitale.
c                   l'indice cont numerote les expressions du type
c                   de l'orbitale en question. ces expressions etant
c                   formees de j (j=1,itypmax) nombres.
c
c
      i=0
 20   i=i+1
 25   e(i)=e(i)+1
      if(e(i).le.e(i-1))go to 30
      e(i)=0
      i=i-1
      if(i.eq.0) go to 10
      go to 25
 30   if(i.lt.ityp) go to 20
      cont=cont+1
c
c
      do 35 j=1,ityp
      per(cont,j)=e(j)
 35   continue
c
c
c
      go to 25
c
c
 100  continue
c
c
c
c *** formation de l'expression numerique decrivant un couple d'orbita-
c     le donne.
c     ..................................................................
c
c
c     noina : numero de la premiere expression numerique de l'orbitale
c             a (ita).
c     noinb : numero de la premiere expression numerique de l'orbitale
c             b (itb).
c
c
      ij=0
      do 400 ii=0,itypmax
      do 300 jj=0,ii
      ij=ij+1
      itypab(1)=ii
      itypab(2)=jj
      do 120 i=1,2
      noin=0
      do 110 j=0,(itypab(i)-1)
      noin=noin+nf(j)
 110  continue
      noin=noin+1
      if(i.eq.1)then
            noina=noin
                else
            noinb=noin
      endif
 120  continue
c
c
c
c     niexp(k) : les 2*itypmax nombres de la ieme integrale sont mis
c                dans le tableau niexp(2*itypmax).
c       itabexp(i,n) : tableau contenant les expressions numeriques qui
c                    decrivent un couple d'orbitale (i=1,nint integra-
c                    les). chaque couple etant defini par n=1,2*itypmax
c                    nombres.
c
c
c
c     ..................................................................
c *** notons que l'ordre des orbitales "d" et "f" qui sont generees dans
c     la subroutine express n'est pas le meme que dans "hondo". le ta-
c     bleau tabor rectifie cet ordre.
c
c     la subroutine express genere les "d" et les "f" respectivement
c     comme suit :
c     xx,yx,yy,zx,zy,zz et xxx,yxx,yyx,yyy,zxx,zyx,zyy,zzy,zzz.
c
c     dans hondo, l'ordre des "d" et des "f" est xx,yy,zz,xy,xz,yz et
c     xxx,yyy,zzz,xxy,xxz,yyx,yyz,zzx,zzy,xyz respectivement.
c     ..................................................................
c
c
c
      i=0
      do 250 j=noina,noina+nf(ii)-1
      do 210 k=1,itypmax
      niexp(k)=per(tabor(j),k)
 210  continue
      do 240 l=noinb,noinb+nf(jj)-1
      do 220 m=1,itypmax
      niexp(itypmax+m)=per(tabor(l),m)
 220  continue
      i=i+1
      do 230 n=1,2*itypmax
      itabexp(ij,i,n)=niexp(n)
 230  continue
 240  continue
 250  continue
 300  continue
 400  continue
c
c
c
      return
      end
