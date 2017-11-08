!
!     DB (06/2014) : 
!       Appels à 'giveis' remplacés par des appels à 'diagonaliser'
!

      subroutine extr(ndet,nv,v,vr,e,qictot,ityp,iiii)
      implicit real*8(a-h,o-p,r-z),logical*4(q)
      character*40 typener
      logical*1 qion
      include 'bd.prm'
c calcul d'hamiltonien effectif dans la base de reference
      integer*4 ne,trpa
      dimension sre(metz,metz),srd(metz,metz),ex(metz),exd(metz)
      dimension z(metz*metz),zd(metz*metz),h0(metz,metz),x(metz,metz)
      dimension h(metz,metz),hd(metz,metz),exbr(metz*metz)	
      dimension exvp(metz*metz),extot(metz*metz),et(metz)
      dimension zucu(metz*metz)
      dimension v(*),vr(metz*2000),e(metz),w(metz*metz)   
c     logical to skip dav. correction if problem
      qfail=.false. 
      
      ndim=20
      write(6,*)   ' extra',nv,ndet,ityp,iii
      qw=.not.((ityp.ne.0).and.(iiii.gt.1))
c
c   reinitialisation e,v
       do i=1,nv
	 et(i)=e(i)
       enddo
       rewind 4
       do i=1,metz
         do j=1,metz
           h0(i,j)=0.d0
           sre(i,j)=0.d0
         enddo
       enddo  
       write(6,*)   '  lecture f04'
c     lecture des donnees
c lecture des etats de references

      read(4,end=99,err=99) nomr,noca,nocb,nvr,ndetr,isz,qion,ii,ik,
     *(ne,j=1,2*ndetr),(trpa,j=1,2*ii),
     *(vr(j),j=1,ik),(h0(j,j),j=1,nvr)          
     
       write(6,*) nvr,ndetr
       write(6,*) 'nomr = ',nomr
      nvbr=nvr*(nvr+1)/2
      read(4)(exvp(i),i=1,nvbr),
     *escf,ityper,ntrsy,(trpa,i=1,nvr*ntrsy)
      if (ityp.ne.0)  ityper = iiii
      do i=1,nvbr
	exvp(i)=1.0
      enddo
c
c si ic totale pas besoin d hamiltonien effectif
c
      if(qictot) then
      write(6,*)
      write(6,*) ' escf ',escf
cdudu
      typener='energie de polarisation des coeurs'
      call rdener(typener,epolnuc,1)
      write(6,*) ' energie de polarisation des coeurs ',epolnuc
cdudu
      write(6,*)
      write(6,*) '    var          var + epolnuc'
c
c     tri
c
      do i=1,nvr
         do j=i+1,nvr
           if(e(i).gt.e(j)) then
                 ej=e(j)
                 e(j)=e(i)
                 e(i)=ej
            endif
          enddo
        enddo    

      do  i=1,nvr
         write(6,985) e(i)+escf,e(i)+escf+epolnuc
      enddo
      stop
      endif
      
      write(6,*) ' extr ityper ',ityper
      if (ityper.eq.1) read(4) (exbr(i),i=1,nvbr)
      if (ityper.eq.2) read(4) (xx,i=1,nvbr),(exbr(i),i=1,nvbr)
      if (ityper.eq.3) read(4) (xx,i=1,2*nvbr),(exbr(i),i=1,nvbr)
      if (ityper.eq.1) read(4) (exvp(i),i=1,nvbr)
      if (ityper.eq.2) read(4) (xx,i=1,nvbr),(exvp(i),i=1,nvbr)
      if (ityper.eq.3) read(4) (xx,i=1,2*nvbr),(exvp(i),i=1,nvbr) 
      nelec=2*noca
      if(qion) nelec=nelec-1
      xelec=nelec
      scale=(xelec-2.d0)/xelec
      if(nelec.le.2) scale=0.d0

c
c
c    lecture des etats a projeter
c
      iv=ndet*nv
      if(qw) write(6,255) escf
255   format(/1x,' energie de la reference couche fermee ',f12.6,/)
c
c  recouvrement entre etats et etats de reference
c  sre: recouvrement moyen/espace de reference cipsi
c       projecteur=somme /psicip><psicip/
c       =>heff
c  srd: recouvrement moyen/espage generateur s
c       projecteur=somme /dets><dets/
c       =>extrapolation q davidson
c
      ndett=min0(ndet,ndetr)
      if(qw) then
      write(6,*)
      write(6,*) ' les references de heff sont les etats de cipsi'
      write(6,*) ' nombre de determinants de reference: ',ndetr
      write(6,*) ' nombre d etats de reference: ',nvr
      write(6,*) ' nombre d etats adiabatiques: ',nv
      endif
       do 825 j=1,nv
       mj=(j-1)*ndet
          am=0.d0
           ad=1.d0
       do 824 i=1,nvr
       mi=(i-1)*ndetr
       mii=(i-1)*ndet
       a=0.d0
         b=0.d0
       do 830 k=1,ndett
         b=b+v(mii+k)*v(mj+k)
830    a=a+vr(mi+k)*v(mj+k)
         if(abs(a).gt.am)then
           am=abs(a)
            ad=sign(ad,a)
           jj=j
           endif
          srd(i,j)=b
824    sre(i,j)=a
        do823i=1,nvr
         sre(i,jj)=sre(i,jj)*ad
823       continue
          do822k=1,ndet
         mj=(jj-1)*ndet
822      v(mj+k)=ad*v(mj+k)
825         continue
      if(qw) then
      write(6,840)
840   format(/,' recouvrement entre etats et references: < r / psi >')
      kk=0
      do 851 l=1,nvr,13
      kk=kk+1
      do 850 i=1,nvr
850   write(6,596)i,(sre(i,j),j=(kk-1)*13+1,min0(nvr,kk*13))
851   continue
596   format(1x,i3,13f9.5)
      endif
c
c     coefficient  davidson
      if(qw) then
      write(6,*)
      write(6,*) ' coefficient de correction davidson '
      endif
      do  i=1,nvr
      ex(i)=1.d0/dsqrt(srd(i,i))
      if(qw)   write(6,*)ex(i)
      enddo

c
       call low(sre,metz)
      nvbb=nv
c
c    recouvrement entre fonctions modeles
c
      ij=0
      do 880 j=1,nv
      do 880 i=1,j
      ij=ij+1
      xx=0.d0
      do 882 k=1,nvr
882   xx=xx+sre(k,i)*sre(k,j)
880   z(ij)=xx
c
c   fonctions modeles orthonormees de decloizeaux (s-1/2)
c
      if(nvr.eq.1) then
      exd(1)=z(1)
      hd(1,1)=1.d0
      else
      !call giveis(nvr,nvr,nvr,z,h,exd,zd,ierr)
      call diagonaliser( nvr, z, exd, zd,nvr )
      if(qw) then
      write(6,*) ' valeurs propres de la matrice de recouvrement'
      kk=0
      do 2851 l=1,nvr,13
      kk=kk+1
      do 2850 i=1,nvr
2850   write(6,1596) (exd(j),j=(kk-1)*13+1,min0(nvr,kk*13))
2851   continue
      write(6,*) ' '
      endif
      kij=0
      do 1884 i=1,nvr
      do 1884 j=1,nvr
      kij=kij+1
1884  hd(j,i)=zd(kij)
      endif
      do 884 i=1,nvr
884   exd(i)=1.d0/dsqrt(exd(i))
      ij=0
      do 886 i=1,nvr
      do 886 j=1,i
      ij=ij+1
      xx=0.d0
      do 888 l=1,nvr
888   xx=xx+hd(i,l)*hd(j,l)*exd(l)
886   z(ij)=xx
      do 890 i=1,nvr
      do 890 j=1,nvr
      xx=0.d0
      do 892 k=1,nvr
      if(i.ge.k) then
      ik=i*(i-1)/2+k
      else
      ik=k*(k-1)/2+i
      endif
892   xx=xx+z(ik)*sre(j,k)
890   hd(j,i)=xx
c
c    heff des cloiseaux
c
      kl=0
      do 702 k=1,nvr
      do 702 l=1,k
      kl=kl+1
      a=0.d0
      do 703 i=1,nvr
703   a=a+hd(i,l)*e(i)*hd(i,k)
      z(kl)=a
702   h(k,l)=a
      if(qw) then
      write(6,706)
706   format(/,1x,'heff   des cloizeaux ')
      do 704 i=1,nvr
704   write(6,595) ((h(i,j)),j=1,i-1),h(i,i)+escf
      write(6,*)
      endif
c
c     verification
c
c
c
      kl=0
      do 1702 k=1,nvr
      do 1702 l=1,k
      kl=kl+1
      a=0.d0
      do 1703 i=1,nvr
1703  a=a+hd(i,l)*hd(i,k)
1702  zucu(kl)=a
      if(qw) then
      write(6,*) ' '
      write(6,*) ' ucroix u:'
      kk=0
      do 1851 l=1,nvr,13
      kk=kk+1
      do 1850 i=1,nvr
1850   write(6,1596) (zucu(j),j=(kk-1)*13+1,min0(nvr,kk*13))
1851   continue
1596   format(1x,13f9.5)
       write(6,*) ' '
       write(6,*) ' '
       endif
c
      !call giveis(nvr,nvr,nvr,z,sre,exd,srd,ierr)
      call diagonaliser( nvr, z, exd, srd,nvr )
      do 7705 i=1,nvr
      if (dabs(e(i)-exd(i)).gt.1d-6) then
      write(6,*) ' attention erreur erreur erreur -----  $*#@%&!??'
      write(6,*) ' heff ne reproduit pas la solution exacte  @$!#?'
      write(6,*)' v. p. de heff',(exd(ki),ki=1,nvr)
      write(6,*)' sol. exacte  ',(e(ki),ki=1,nvr)
      qfail=.true.
      goto 7706
      endif
7705  continue
      write(6,*) ' heff reproduit la solution exacte'
7706  continue
c
c  davidson
c
        do460i=1,nv
          do460j=1,i
460    hd(i,j)=(h(i,j)-h0(i,j))*(ex(i)*ex(j)-1.d0)
c
c   interpolation lineaire sur le reste de perturbation
c
      write(6,*)
      write(6,*) ' reste de la perturbation (ityper = ',ityper,')'
      write(6,*) ' '
      if (ityper .eq. 1) write(6,*) ' CALCUL MOLLER PLESSET
     * BARYCENTRIQUE'
      if (ityper .eq. 2) write(6,*) ' CALCUL EPSTEIN NESBET
     * VALEUR PROPRE'
      if (ityper .eq. 3) write(6,*) ' CALCUL EPSTEIN NESBET
     * BARYCENTRIQUE'
      write(6,*)  ' ' 

      ij=0
      do 887 i=1,nvr
      do 889 j=1,i
      ij=ij+1
      vr(ij)=exbr(ij)
      x(i,j)=exbr(ij)-h0(i,j)
889   if(.not.qictot.or.j.eq.i) x(i,j)=x(i,j)-exvp(ij)
      x(i,i)=x(i,i)-escf
      write(6,595) (x(i,j),j=1,i)
887   continue

c
c*************************************
c     hamiltoniens effectifs divers
c**************************************
c
c   resultats precedents
c
      call diagonaliser( nvr, exbr, ex, srd,nvr )
      write(6,*)
      write(6,*) ' resultats de cipsi Partition S'
      write(6,*)
                
         write(6,*)'  ordre zero     ordre deux  '
      do i=1,nvr
      write(6,985) h0(i,i)+escf,ex(i)
      enddo

      ij=0
      do i=1,nvr
      do j=1,i-1
      ij=ij+1
      extot(ij)=x(i,j)
      enddo
      ij=ij+1
      extot(ij)=e(i)+escf+x(i,i)
      enddo
      call diagonaliser( nvr, extot, ex, srd,nvr )
      write(6,*)
      write(6,*) ' resultats de cipsi (diagonalisation de heff)'
      write(6,*)
                
      write(6,*)'  ordre zero     ordre deux  '
      do i=1,nvr
      write(6,985) e(i)+escf,ex(i)
      enddo
      
c     RV2016 : return if failure 
      if(qfail) then
        do i=1,nvr
        e(i)=ex(i)
        enddo
        write(6,*)
        write(6,*) ' escf ',escf
        typener='energie de polarisation des coeurs'
        call rdener(typener,epolnuc,1)
      write(6,*) ' energie de polarisation des coeurs ',epolnuc
      write(6,*)
      write(6,*) '    var          var + epolnuc'
      do i=1,nvr
         do j=i+1,nvr
           if(e(i).gt.e(j)) then
                 ej=e(j)
                 e(j)=e(i)
                 e(i)=ej
            endif
          enddo
        enddo    
      do  i=1,nvr
         write(6,985) e(i),e(i)+epolnuc
      enddo
      return
      endif
      
c  resultats sans correction  davidson
c
      nbr1=nvbr
      nbr2=nbr1+nvbr
      ij=0
      do 952 i=1,nvr
      do 953 j=1,i
      ij=ij+1
      z(ij)=h(i,j)+x(i,j)
      if(dabs(exvp(ij)).gt.1.d-10) then
      zd(ij)=h(i,j)+x(i,j)*(h(i,j)-h0(i,j))/exvp(ij)
      else
      zd(ij)=z(ij)
      endif
      vr(nbr1+ij)=z(ij)
      vr(nbr2+ij)=zd(ij)
953   continue
      z(ij)=z(ij)+escf
      zd(ij)=zd(ij)+escf
      vr(nbr1+ij)=z(ij)
      vr(nbr2+ij)=zd(ij)
952   continue
      write(6,*) ' heff var+reste'
      ndec=0
      do 1704 i=1,nvr
      write(6,595) (z(ndec+j),j=1,i)
      ndec=ndec+i
1704  continue
      write(6,*)
      write(6,*)
cx      write(6,*) ' heff var+reste extra'
cx     ndec=0
cx      do 1705 i=1,nvr
cx      write(6,595) (zd(ndec+j),j=1,i)
cx      ndec=ndec+i
cx 1705 continue
cx      write(6,*)
       !call giveis(nvr,nvr,nvr,z,sre,ex,srd,ierr)
       !call giveis(nvr,nvr,nvr,zd,sre,exd,srd,ierr)
       call diagonaliser( nvr, z, ex, srd,nvr )
       call diagonaliser( nvr, zd, exd, srd,nvr )
       write(6,*)
       write(6,*) '***************************************'
       write(6,*) '  energies apres ic des moyens  (s+m)  '
       write(6,*) '***************************************'
       write(6,*)
      write(6,*) '   var        var+reste      var+reste extra'
      ij=0
      do 954 i=1,nvr
      ij=ij+1
954   write(6,985) e(i)+escf,ex(i),exd(i)
c
c  davidson independant du nombre d electrons
c
      nbr1=nbr2+nvbr
      nbr2=nbr1+nvbr
      write(6,*)
      write(6,*) ' correction davidson  usuelle q '
      ij=0
      do 956 i=1,nvr
      do 957 j=1,i
      ij=ij+1
      w(ij)=h(i,j)+hd(i,j)
      z(ij)=w(ij)+x(i,j)
      if(dabs(exvp(ij)).gt.1.d-10) then
      zd(ij)=w(ij)+x(i,j)*(h(i,j)+hd(i,j)-h0(i,j))/exvp(ij)
      else
      zd(ij)=z(ij)
      endif
      vr(nbr1+ij)=z(ij)
      vr(nbr2+ij)=zd(ij)
957   continue
      w(ij)=w(ij)+escf
      z(ij)=z(ij)+escf
      zd(ij)=zd(ij)+escf
      vr(nbr1+ij)=z(ij)
      vr(nbr2+ij)=zd(ij)
956   continue
!        call giveis(nvr,nvr,nvr,w,sre,ex,srd,ierr)
!        call giveis(nvr,nvr,nvr,z,sre,exd,srd,ierr)
!        call giveis(nvr,nvr,nvr,zd,sre,e,srd,ierr)
       call diagonaliser( nvr, w, ex, srd,nvr )
       call diagonaliser( nvr, z, exd, srd,nvr )
       call diagonaliser( nvr, zd, e, srd,nvr )
      write(6,*)
      write(6,*) '   var+q      var+q+reste     var+q+reste extra'
      do  i=1,nvr
         write(6,985) ex(i),exd(i),e(i)
      enddo
c
c  sdtq renormalise au nombre d electrons
c
      write(6,969) nelec
969   format(/,' correction q normalisee au nombre d electrons ',i4)
      nbr1=nbr2+nvbr
      nbr2=nbr1+nvbr
      ij=0
      do 966 i=1,nvr
      do 967 j=1,i
      ij=ij+1
      w(ij)=h(i,j)+scale*hd(i,j)
      z(ij)=w(ij)+x(i,j)
      if(dabs(exvp(ij)).gt.1.d-7) then
      zd(ij)=w(ij)+x(i,j)*(h(i,j)+scale*hd(i,j)-h0(i,j))/exvp(ij)
      else
      zd(ij)=z(ij)
      endif
      vr(nbr1+ij)=z(ij)
      vr(nbr2+ij)=zd(ij)
967   continue
      w(ij)=w(ij)+escf
      z(ij)=z(ij)+escf
      zd(ij)=zd(ij)+escf
      vr(nbr1+ij)=z(ij)
      vr(nbr2+ij)=zd(ij)
966   continue
!        call giveis(nvr,nvr,nvr,w,sre,ex,srd,ierr)
!        call giveis(nvr,nvr,nvr,z,sre,exd,srd,ierr)
!        call giveis(nvr,nvr,nvr,zd,sre,e,srd,ierr)
       call diagonaliser( nvr, w, ex, srd,nvr )
       call diagonaliser( nvr, z, exd, srd,nvr )
       call diagonaliser( nvr, zd, e, srd,nvr )
      write(6,*)
      write(6,*) '   var+q      var+q+reste     var+q+reste extra'
      nvbr=nvr*(nvr+1)/2
      do 968 i=1,nvr
968   write(6,985) ex(i),exd(i),e(i)
      rewind 4
      write(6,*) ityper
      do i=1,4
         read(4)
      enddo
      nheff=7
      write(4) (vr(ii),ii=1,nheff*nvbr)
      write(6,*) 'ecriture sur f04'
      write(6,977)
977   format(//,1x,100('*'))
985   format(1x,3(f20.14,2x))
595   format(12f11.6)

c     ecriture final      
      write(6,*)
      write(6,*) ' escf ',escf
cdudu
      typener='energie de polarisation des coeurs'
      call rdener(typener,epolnuc,1)
      write(6,*) ' energie de polarisation des coeurs ',epolnuc
cdudu
      write(6,*)
      write(6,*) '    var          var + epolnuc'
c
c     tri
c
      do i=1,nvr
         do j=i+1,nvr
           if(e(i).gt.e(j)) then
                 ej=e(j)
                 e(j)=e(i)
                 e(i)=ej
            endif
          enddo
        enddo    

      do  i=1,nvr
         write(6,985) e(i),e(i)+epolnuc
      enddo

	 do i=1,nv
	   e(i)=et(i)
         enddo
         return
99      write(6,*) 'error reading file f04'    
        stop
            end
