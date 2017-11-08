      subroutine trvm(xv,xm,idv,idm,jdm,iamin,iamax,jamin,jamax,sym)
      implicit real*8 (a-h,o-z)
      logical sym
      dimension xm(idm,jdm)
      dimension xv(idv)
c.....transfert du vecteur xv dans la matrice xm
c..... xm est stockee en colonnes
c.....idv  dimension de xv
c.....idm,jdm dimensions de xm
c.....iamin,iamax bornes de colonne de xm
c.....jamin,jamax bornes de lignes de xm
c.....sym=.true. la matrice xm est symetrique et xv ne contient que
c.....la demi-matrice inferieure
c.....sym=.false. la matrice xm n'est pas symetrique
      if(sym.and.(iamin.ne.jamin.or.iamax.ne.jamax))then
	write(6,*)' erreur dans trvm ,dans le cas symetrique les bornes
     < de la matrice doivent etre identiques'
	stop
      end if
      if(sym) then
        ij=0
        do i=iamin,iamax
          do j=jamin,i 
	    ij=ij+1
	    xm(i,j)=xv(ij)
	    xm(j,i)=xv(ij)
          end do
        end do
      else
        ij=0
        do j=jamin,jamax
          do i=iamin,iamax
	    ij=ij+1
	    xm(i,j)=xv(ij)
          end do
        end do
      end if
      return
      end
