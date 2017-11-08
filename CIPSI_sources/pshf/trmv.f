      subroutine trmv(xm,xv,idm,jdm,idv,iamin,iamax,jamin,jamax,sym)
      implicit real*8 (a-h,o-z)
      logical sym
      dimension xm(idm,jdm)
      dimension xv(idv)
c.....transfert de la matrice xm dans le vecteur xv
c..... xm est stockee en colonnes
c.....idv  dimension de xv
c.....idm,jdm dimensions de xm
c.....iamin,iamax bornes de colonne de xm
c.....jamin,jamax bornes de lignes de xm
c.....sym=.true. la matrice xm est symetrique et xv ne contient que
c.....la demi-matrice inferieure
c.....sym=.false. la matrice xm n'est pas symetrique
      if(sym.and.(iamin.ne.jamin.or.iamax.ne.jamax))then
	write(6,*)' erreur dans trmv dans le cas symetrique les bornes
     < de la matrice doivent etre identiques'
	stop
      end if
      if(sym) then
        ij=0
        do i=iamin,iamax
          do j=jamin,i 
	    ij=ij+1
	    xv(ij)=xm(i,j)
          end do
        end do
      else
        ij=0
        do i=iamin,iamax
          do j=jamin,jamax
	    ij=ij+1
	    xv(ij)=xm(i,j)
          end do
        end do
      end if
      return
      end
