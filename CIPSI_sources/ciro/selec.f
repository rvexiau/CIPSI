      subroutine selec(yuni,mspin,typ_sym1,lab1,metat1,
     *                            typ_sym2,lab2,metat2,
     *                            ietat,jetat,ncou)
      implicit real*8 (a-h,o-x,z),logical*1 (y)
      include 'pshf.prm'
      character*40 typener,typ_sym1,typ_sym2
      character*3 lab1,lab2,sym
      dimension ietat(*),jetat(*),istate(metz),sym(metz)
      dimension num1(metz),num2(metz)
      rewind 33
      nenr=0
    1 read(33,end=9000)typener
      nenr=nenr+1
      if(typener.ne.typ_sym1) go to 1
      rewind(33)
      do i=1,nenr-1
	read(33)
      end do
      read(33)typener,(istate(i),sym(i),i=1,metat1)
      ns1=0
      do i=1,metat1
	if(istate(i).eq.mspin.and.sym(i).eq.lab1)then
	  ns1=ns1+1
	  num1(ns1)=i
	end if
      end do
      if(yuni)then
	ns2=ns1
	do i=1,ns2
	num2(i)=num1(i)
	end do
	else
      rewind 33
      nenr=0
    2 read(33,end=9000)typener
      nenr=nenr+1
      if(typener.ne.typ_sym2) go to 2
      rewind(33)
      do i=1,nenr-1
	read(33)
      end do
      read(33)typener,(istate(i),sym(i),i=1,metat2)
      ns2=0
      do i=1,metat2
	if(istate(i).eq.mspin.and.sym(i).eq.lab2)then
	  ns2=ns2+1
	  num2(ns2)=i
	end if
      end do
      end if
      if(ns1.eq.0.or.ns2.eq.0) then    
	write(6,*) ' pas d''etat correct de symetrie ',mspin, ' sur la
     * file 33'
      stop
      end if
      write(6,*) ns1,' etats',mspin,lab1,' et',ns2,' etats',mspin,lab2, 
     * '  etats selectiones sur la file 33'
      do i=1,ns1
	do j=1,ns2
	  ncou=ncou+1
	  if(ncou.gt.ncouz) then
	    write(6,*) ' ncouz trop petit dans ciro,stop'
	    stop
          end if
	  ietat(ncou)=num1(i)
	  jetat(ncou)=num2(j)
        end do
      end do
      write(6,*)' selection automatique des etats multip=  ',mspin,
     * 'sym1  ',lab1, '  sym2',lab2
      write(6,*)' nombre de couples', ncou
      return
 9000 write(6,*)' erreur sur file ener (33) pas de valeurs pour le label
     * ', typ_sym1,typ_sym2
      stop
      end
