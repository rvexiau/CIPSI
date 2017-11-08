      subroutine matbc(a,ida,jda,iamin,iamax,jamin,jamax,
     <                 b,idb,jdb,ibmin,ibmax,jbmin,jbmax,
     <                 c,idc,jdc,icmin,icmax,jcmin,jcmax,
     <                 x,idx,jdx,ixmin,ixmax,jxmin,jxmax,
     <                 temp)
c.....                     t
c.....calcul du produit x=a*b*c
c.....ida,jda dimensions de a
c.....idb,jdb dimensions de b
c.....idc,jdc dimensions de c
c.....idx,jdx dimensions de x
c.....iamin,iamax bornes de l'indice de colonne de a
c.....jamin,jamax bornes de l'indice de ligne de a
c.....ibmin,ibmax bornes de l'indice de colonne de b
c.....jbmin,jbmax bornes de l'indice de ligne de b
c.....icmin,icmax bornes de l'indice de colonne de c
c.....jcmin,jcmax bornes de l'indice de ligne de c
c.....ixmin,ixmax bornes de l'indice de colonne de x
c.....jxmin,jxmax bornes de l'indice de ligne de x
      implicit real*8 (a-h,o-z)
      dimension a(ida,jda),b(idb,jdb),c(idc,jdc),x(idx,jdx)
      dimension temp(idb)
c.....verification des bornes
      if(jamax-jamin.ne.ixmax-ixmin) then
	write(6,*) ' erreur matbc jamin,jamax,ixmin,ixmax',
     <                            jamin,jamax,ixmin,ixmax
	stop
      end if
      if(jcmax-jcmin.ne.jxmax-jxmin) then
	write(6,*) ' erreur matbc jcmin,jcmax,jxmin,jxmax',
     <                            jcmin,jcmax,jxmin,jxmax
	stop
      end if
      if(iamax-iamin.ne.ibmax-ibmin) then
	write(6,*) ' erreur matbc iamin,iamax,ibmin,ibmax',
     <                            iamin,iamax,ibmin,ibmax
	stop
      end if
      if(icmax-icmin.ne.jbmax-jbmin) then
	write(6,*) ' erreur matbc icmin,icmax,jbmin,jbmax',
     <                            icmin,icmax,jbmin,jbmax
	stop
      end if
c
c
      jxx=jxmin
      do jc=jcmin,jcmax
	do ib=ibmin,ibmax
	  dum=0.d0
	  icx=icmin
	  do jb=jbmin,jbmax
	    dum=dum+b(ib,jb)*c(icx,jc)
	    icx=icx+1
          end do
	  temp(ib)=dum
        end do
c
	ixx=ixmin
	do ja=jamin,jamax
	  dum=0.d0
	  ibx=ibmin
	  do ia=iamin,iamax
	    dum=dum+a(ia,ja)*temp(ibx)
	    ibx=ibx+1
          end do
	  x(ixx,jxx)=dum
	  ixx=ixx+1
        end do
	jxx=jxx+1
      end do
      return
      end
