      program cip                                                    
      implicit real*8(a-h,o-x,z),logical*1 (y)                          
      include 'pshf.prm'

      character *40 typener
      character*20 title
      character*26 timeday
      character*4 ibla2
      character*8 group                                                 
      integer*4 ne,nd,trou,part,iorb,iwspin,itsym,its,isytr,nesp,nespo  
      integer*4 itm                                       
      dimension nps(8),nts(8)                                           
      dimension estk(300)
      dimension empdet(ndetz)                                           
      dimension vect(metz*metz)
      character*7 tabsym                                                
      character*3 typsym(metz),sym 
      dimension ispin_etat(metz)
      dimension fmul(metz),tabsym(10)                                   
      common/recla/iord(ndetz),mspin(metz)                              
      common/grou/group                                                 
      common/hist/nom(10),nam(10,metz),aval(10),tabmp(10,metz),         
     1tau,semp(metz),sevp(metz),tabvp(10,metz),nclass,ymoyen            
      common/det/my,ny,idet(400),nexst(100),vst(100),long(8)            
      common/spo/fdiag(2*doa),hdiag(ndetz),tocom(2*doa),coeff,      
     1ndef,ispin(2*doa),iwspin(2*doa),numero(metz),          
     2ne(ndetz+1),nd(ndetz+1),trou(ntpz),part(ntpz),iorb(2*doa),      
     3itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),yocs(2*doa),       
     4ygen(ndetz),yocd(2*doa,ndetz),yspin                                    
      common/degene/nesp(ndetz),nespo(ndetz),isytr(ndetz,nsymz),        
     1yef(ndetz,nsymz),ysaut(ndetz)                                     
      common/nature/itm(metz,nsymz),ityper                              
      common/tome/taum(metz),yselm,ysele,ydp                            
      common/nom33/title
      common/ij1/ij0(doa,doa),ntttt(nsymz)      
      common/cipsi/escf,test,teseff,tdet,twdet,ndiag,norb,mnorb,
     1             noca,nocb,noc1,noc2,igela,igelb,ncf,nsym,isym,isz,
     2             ntrsy,metat,lecdet,maxdet,yprt,yion,ypertu,ybrd,yen,
     3             ywvi,ywvf,ytoul,ystkic,ywdeg
      common/ic/c(metz*ndetz),e(ndetz),eb(metz),e2b(metz),emp(metz),    
     1teste(metz),e2mp(metz),dpivot,ietat,ishif,e2vp(metz),             
     2fmpb(2*doa),hefmp(nhefz),hefvp(nhefz),hefb(nhefz),              
     3psimp(metz),psivp(metz),hzero(nhefz)                              
      common/giv/nblock(ndetz),nstart(ndetz)                            
      data ibla1,ibla2/0,'=>'/
      call fdate(timeday)
       write(6,9999)timeday
 9999 format(/,80(1h*),/,8(10h  cipsi   ),/,  10x,a25,/,
     * 80(1h*))
      call openf                                                        
      do 205 i=1,2*doa                                                
205   tocom(i)=0.d0                                                     
c                                                                       
      call deter                                                        
      nesp(ncf)=0                                                       
c                                                                       
      noc1=max0(noca,nocb-norb)                                         
      noc2=min0(noca,nocb-norb)                                         
      call ijkf
      write(*,*) 'lecture des integrales '
      call reijkl(norb,nsym,ydp)                                        
      ne1=ne(1)                                                         
      nd1=nd(1)                                                         
      call initsm(tabsym,nsym)                                          
      write(6,8001) group                                               
8001  format(1x,'groupe de symetrie: ',a8)                              
      write(6,*)    
      typener='energie du determinant de reference'
      call rdener(typener,escf,1)      
      if(ityper.ge.0) call ic001                                                        
      do 8 i=1,mnorb                                                    
8      yoc(i)=.false.                                                   
c                                                                       
c   calcul de l operateur de fock                                       
      do 760 i=1,norb                                                   
      if(tocom(i).gt.1.d0) go to 755                                    
      tocom(i+norb)=0.d0                                                
      go to 760                                                         
755    tocom(i+norb)=tocom(i)-1.d0                                      
      tocom(i)=1.d0                                                     
760   continue                                                          
      do 720 i=1,noca                                                   
      tocom(i)=tocom(i)-1.d0                                            
720    tocom(i+norb)=tocom(i+norb)-1.d0                                 
      call faifoc                                                       
      do 790  i=1,norb                                                  
      fmpb(i)=(fmpb(i)+fmpb(i+norb))*0.5d0                              
790   fmpb(i+norb)=fmpb(i)                                              
      write(6,4600) (fmpb(l),l=1,norb)                                  
4600  format(1x,'energies monoelectroniques',//,10(1x,f8.3))            
c                                                                       
c   energies d ordre zero des etats                                     
c                                                                       
c   energies mp des determinants                                        
c                                                                       
      write(6,1105) escf                                                
1105  format(1x,'energie de reference:',f15.6,//)                       
c
      if(metat.eq.1) ybrd=.false.                                       
      do 25 i=1,ncf                                                     
25    empdet(i)=hmp(i)                                                  
      jn=-ncf                                                           
      do 21 m=1,metat                                                   
      e2vp(m)=0.d0                                                      
      e2mp(m)=0.d0                                                      
      e2b(m)=0.d0                                                       
      emp(m)=0.d0                                                       
      t=0                                                               
      jn=jn+ncf                                                         
      do 20 i=1,ncf                                                     
      ij=jn+i                                                           
      km=ij                                                             
      if(ysaut(i).or.t.gt.dabs(c(km)))  go to 18                        
      t=dabs(c(km))                                                     
      nrk=i                                                             
18    continue                                                          
c                                                                       
c  moller plesset                                                       
c                                                                       
      emp(m)=emp(m)+c(ij)*c(ij)*empdet(i)                               
20    eb(m)=eb(m)+(c(ij)*c(ij)*hdiag(i))                                
c                                                                       
c  les seuils reels de selection des determinants pour chaque etat      
c  est ajuste sur le plus grand coefficient de la fonction              
c  multiconf de cet etat                                                
c                                                                       
      teste(m)=test                                                     
      taum(m)=1.d0                                                      
      if(yselm)then                                                     
        teste(m)=t*test                                                 
       if(ityper.ge.0) then
        taum(m)=1.d0/t                                                  
       endif
      endif                                                             
21    continue                                                          
c                                                                       
c     ecriture des caracteristiques d'ordre zero de la fonction sur la  
c     file 4                                                            
      jj=ncf*metat                                                       
      ii=nd(ncf)+ne(ncf)+1                                              
      rewind 4                                                          
      e0=0.d0
      if(ityper.eq.0) e0=escf
      write(4) norb,noca,nocb,metat,ncf,isz,yion,ii,jj,(ne(j),nd(j),j=1  
     * ,ncf),                                                           
     y (trou(j),part(j),j=1,ii),(c(j),j=1,jj),(e(j)+e0,j=1,metat),          
     y (emp(j),j=1,metat),(eb(j),j=1,metat),(numero(i),i=1,metat),      
     y (empdet(i),i=1,ncf),(fmpb(j),j=1,norb),(tocom(j),j=1,norb),      
     y (nblock(i),i=1,ndetz),(nstart(i),i=1,ndetz),(iorb(i),i=1,2*norb)    
     y ,(ispin(i),i=1,2*norb),((yocd(i,j),i=1,2*norb),j=1,ncf)    
  115 continue                                                          
      write(6,1030)                                                     
      write(6,1030)                                                     
      write(6,1010)                                                     
1010  format(1x,'seuil effectif de selection et energies d ordre zero') 
      write(6,1030)                                                     
      if(.not.yen) then                                                 
      write(6,1005)                                                     
1005  format(1x,'etat',5x,'test',8x,'tau',8x,'e0-mpb')                  
      else                                                              
      write(6,1020)                                                     
1020  format(1x,'etat',5x,'test',8x,'tau',8x,'e0-mpb',                  
     * 6x,'e0-envp',7x,'e0-enb')                                        
      endif                                                             
1030  format(/)                                                         
      do 11 i=1,metat                                                   
      if(.not.yen) then                                                 
      write(6,1035) numero(i),teste(i),tau/taum(i),emp(i)               
1035  format(2x,i2,5(4x,f10.6))                                          
      else                                                              
      write(6,1035) numero(i),teste(i),tau/taum(i),emp(i),e(i),eb(i)    
      endif                                                             
11    continue                                                          
      write(6,1030)                                                     
333   continue                                                          
      write(6,*) 'avant brdneee'
      call brdnee                                                       
      do 50 m=1,metat                                                   
      e2mp(m)=0.                                                        
      psimp(m)=0.d0                                                     
      psivp(m)=0.d0                                                     
      e2vp(m)=0.0                                                       
50    e2b(m)=0.0                                                        
      do 52 i=1,nhefz                                                   
      hefmp(i)=0.d0                                                     
      hefvp(i)=0.d0                                                     
52    hefb(i)=0.d0                                                      
      ietat =0                                                          
      dpivot =0.0                                                       
      do 40 i=1,mnorb                                                   
40     yoc(i)=.false.                                                   
      nf=ncf+1                                                          
      nd(nf)=nd(ncf)+ne(ncf)+1                                          
      nbr=metat*(metat+1)/2                                             
      
c
c si ityper<0 tri des determinants par energies croissantes
c
      if(ityper.lt.0) then
       do  ii=1,ncf
           empdet(ii)=-hdig(ii)
       enddo
       write(6,*) 'avant shell'
       call shell(empdet,iord,ncf)
       write(6,*) 'apres shell'
        do  kk=1,ncf
         k=iord(kk)
         k1=1+nd(k)
         k2=k1+ne(k)-1
       if(ywvi) write(6,6666) (iorb(trou(i)),iwspin(trou(i)),i=k1,k2),
     * ibla1,ibla2, (iorb(part(i)),iwspin(part(i)),i=k1,k2)
6666  format(1x, 8(i3,a2,' '),/,(1x,8(i3,a2,' ')))
       enddo
      
       write(6,*) 'ecriture sur 4'
       write(4) (hzero(ij),ij=1,nbr),                                 
     *   escf,ityper,ntrsy,((itm(m,l),m=1,metat),l=1,ntrsy)             
         write(4) (hefmp(ii),ii=1,nbr),                                 
     1         (hefvp(ii),ii=1,nbr),                                    
     2         (hefb(ii),ii=1,nbr)   
      endif
      if(.not.ymoyen) go to 5                                           
      do 4 kss=1,ncf                                                    
      ks=iord(kss)                                                      
      necst=ne(ks)                                                      
      necsv=necst                                                       
      if(necsv.eq.0) necsv=1                                            
      nposte=nd(ks)                                                     
      do 3 ls=1,necsv                                                   
      nts(ls)=trou(nposte+ls)                                           
3     nps(ls)=part(nposte+ls)                                           
      call pie(necst,nts,nps,1.d0)                                      
4     continue                                                          
      call fine                                                         
5     continue                                                          
      if(ityper.lt.0) stop
      write(6,*) ' avant pertu yreduc=', yreduc
      if(ypertu)call pertu                                              
      write(6,*) ' apres pertu yreduc=', yreduc
      write(6,*) 'fin pertu'
      if(ypertu)call diexcit                                            
      if(ymoyen) call fine                                              
      write(6,1225) ietat                                               
1225   format('0 ont ete generes:    ',i10,'  etats efficaces')         
      write(6,1370) dpivot                                              
1370  format(1x,'la contribution du plus grand determinant',            
     1' non ecrit est ',f8.6,/)                                         
      write(6,1400)                                                     
1400  format('0norme de la correction  d ordre 1 a la fonction',/)      
      if(yen) then                                                      
      write(6,*)  '  etat          mpb           envp'                  
      else                                                              
      write(6,*)  '  etat          mpb'                                 
      endif                                                             
      do 1410 ii=1,metat                                                
      if(yen) then                                                      
      write(6,1085) numero(ii),psimp(ii),psivp(ii)                      
      else                                                              
      write(6,1085) numero(ii),psimp(ii)                                
      endif                                                             
1085  format(2x,i2,2x,f15.6,2x,f15.6)                                   
1410  continue                                                          
      write(6,1030)                                                     
      write(6,2000)                                                     
      call spasym(metat)                                                
      call mult(fmul,metat)                                             
      write(6,1060)                                                     
1060  format(1x,'contributions a l energie totale')                     
      write(6,2000)                                                     
2000  format(1x,32('*'))                                                
      write(6,1030)                                                     
      if(.not.yen) then                                                 
      write(6,1070)                                                     
1070  format(1x,'etat  symetrie',6x,'  diag',6x,' pertu mpb')           
      else                                                              
      write(6,1065)                                                     
1065  format(1x,'etat  symetrie',6x,'  diag',6x,' pertu mpb',           
     * 6x,' pertu envp',6x,' pertu enb')                                
      endif                                                             
      do 10 m=1,metat                                                   
      ks=int(0.5+2.*cnvmul(fmul(m)))                                    
      call recsym(m,sym)                                                
      if(.not.yen) then                                                 
      write(6,1086) numero(m),ks+1,sym,e(m),e2mp(m)                     
1086  format(2x,i2,4x,i1,a3,4(2x,f14.8))                                
      else                                                              
      write(6,1086) numero(m),ks+1,sym,e(m),e2mp(m),e2vp(m),e2b(m)      
      endif                                                             
  10  continue                                                          
      write(6,1030)                                                     
      write(6,1030)                                                     
      if(.not.yen) then                                                 
      write(6,1110)                                                     
1110  format(1x,'etat  symetrie',4x,' scf+diag',6x,'total mpb')         
      else                                                              
1115  write(6,1120)                                                     
1120  format(1x,'etat  symetrie',4x,' scf+diag',6x,'total mpb',         
     * 7x,'total envp',6x,'total enb')                                  
      endif                                                             
      mm=0                                                              
      do 690 m=1,metat                                                  
      mm=mm+m                                                           
      e(m)=e(m)+escf                                                    
      hefmp(mm)=e(m)+e2mp(m)                                            
      hefvp(mm)=e(m)+e2vp(m)                                            
      hefb(mm)=e(m)+e2b(m)                                              
      ks=int(0.5+2.*cnvmul(fmul(m)))                                    
      call recsym(m,sym)                                                
      typsym(m)=sym
      ispin_etat(m)=ks+1
      if(.not.yen) then                                                 
      write(6,1086) numero(m),ks+1,sym,e(m),hefmp(mm)                   
      else                                                              
      write(6,1086) numero(m),ks+1,sym,e(m),hefmp(mm),                  
     1 hefvp(mm),hefb(mm)                                               
      endif                                                             
690   continue                                                          
      if(ybrd) then                                                     
      write(6,*)                                                        
      write(6,2100)                                                     
      write(6,1300)                                                     
1300  format(1x,'hamiltonien effectif qdpt ')                           
      write(6,2100)                                                     
2100  format(1x,46('*'))                                                
      write(6,1030)                                                     
      write(6,*) '  heff mpb'                                           
      ijfin=0                                                           
      do 506 i=1,metat                                                  
      ijfin=ijfin+i                                                     
      ijdeb=ijfin-i+1                                                   
506   write(6,1350) (hefmp(ii),ii=ijdeb,ijfin)                          
      call jacscf(hefmp,vect,e,metat,0,1.d-12)                  
      do i=1,metat
         e(i)=-e(i)
      enddo
c
c   ecriture vecteurs
c
      write(6,*)
      write(6,*) ' vecteurs qdpt mpb'
      write(6,*)
      do m=1,metat
       write(6,1350) (vect((m-1)*metat+k),k=1,metat)
      enddo
c
c tri par energies croissantes
c
      call shell(e,iord,metat)
      write(6,*)
      write(6,*) 'ordre des valeurs propres'
      write(6,*)
      write(6,*) (iord(m),m=1,metat)
      write(6,*)
      write(6,*) 'energies heff mpb'                                           
      do  m=1,metat                                                   
         mm=iord(m)
         ks=int(0.5+2.*cnvmul(fmul(mm)))                                    
         call recsym(mm,sym)                                                
         write(6,1086) numero(m),ks+1,sym,-e(m)
      enddo
c
c  recombinaison et stockage des vecteurs QDPT de l espace modele 
c  sur les determinants 
c
      do i=1,ncf
      do mm=1,metat
      m=iord(mm)
      hdiag(m)=0.d0
        do k=1,metat
          hdiag(m)=hdiag(m)+vect((m-1)*metat+k)*c((k-1)*ncf+i)
        enddo
       enddo
      do m=1,metat
        c((m-1)*ncf+i)=hdiag(m)
      enddo
      enddo
c
c  impression des vecteurs propres finaux
c
c
      write(59) ncf,metat,(c(i),i=1,ncf*metat),(e(i),i=1,metat)
      if(yen) then                                                      
      write(6,1030)                                                     
      write(6,*) '  heff envp'                                          
      ijfin=0                                                           
      do 507 i=1,metat                                                  
      ijfin=ijfin+i                                                     
      ijdeb=ijfin-i+1                                                   
507   write(6,1350) (hefvp(ii),ii=ijdeb,ijfin)                          
      write(6,1030)                                                     
      write(6,*) '  heff enb'                                           
      ijfin=0                                                           
      do 508 i=1,metat                                                  
      ijfin=ijfin+i                                                     
      ijdeb=ijfin-i+1                                                   
508   write(6,1350) (hefb(ii),ii=ijdeb,ijfin)                           
      endif                                                             
      endif                                                             
      zero=0.d0                                                         
c                                                                       
c   ecriture sur la file 4  des resultats apres perturbation            
c                                                                       
      nheff=1                                                           
	 do 1345 i=1,nbr                                                       
	    hzero(i)=0.d0                                                      
 1345    continue                                                       
	 do 1348 i=1,metat                                                     
	    hzero(i*(i+1)/2)=e(i)-escf                                         
 1348    continue                                                       
         write(4) (hzero(ij),ij=1,nbr),                                 
     *   escf,ityper,ntrsy,((itm(m,l),m=1,metat),l=1,ntrsy)             
         write(4) (hefmp(ii),ii=1,nbr),                                 
     1         (hefvp(ii),ii=1,nbr),                                    
     2         (hefb(ii),ii=1,nbr)                                      
1350  format(1x,10(1x,f11.6))                                           
      typener='symetries CIPSI '//title
      call stktyp(typener,ispin_etat,typsym,metat)
      typener='energies CIPSI var '//title
      call stkener(typener,e,metat)
      typener='energies CIPSI mpb '//title
      call stkener(typener,e2mp,metat)
      typener='energies CIPSI envp'//title
      call stkener(typener,e2vp,metat)
      typener='energies CIPSI enb '//title
      call stkener(typener,e2b,metat)
      write(6,1030)                                                     
      write(6,1030)                                                     
      if(.not.ymoyen) go to 1800                                        
      write(6,1500)                                                     
      write(6,1510)                                                     
1510  format(1x,'histogramme de la perturbation')                       
      write(6,1500)                                                     
1500  format(1x,30('*'))                                                
      write(6,1030)                                                     
      write(6,1038) teseff                                              
1038  format(1x,' seuil de tri teseff= ',f8.6)                          
      if(yen) then                                                      
      write(6,*) ' tri envp'                                            
      else                                                              
      write(6,*) ' tri mpb'                                             
      endif                                                             
      if(ysele) then                                                    
      write(6,*) ' critere: contribution a l energie'                   
      else                                                              
      write(6,*) ' critere: contribution a la fonction'                 
      if(yselm) write(6,*)' critere renormalise  par etat '             
      endif                                                             
      write(6,1030)                                                     
      write(6,1520)                                                     
1520  format(1x,'inf est la borne inf de chaque classe',/,              
     *1x,'sup est la borne sup de chaque classe',/,                     
     *1x,'ntotal est le nombre total de determinants engendres',        
     *' dans chaque classe ',/,                                         
     *1x,'n est le nombre de determinants engendres dans chaque',       
     *' classe pour chaque etat',/,                                     
     *1x,'emp est la contribution mpb  de chaque classe ',              
     *'pour chaque etat',/,                                             
     *1x,'evp est la contribution envp de chaque classe',               
     *' pour chaque etat',//)                                           
      write(6,1530) (ii,ii=1,nclass)                                    
1530  format(1x,'classes:',8x,10(4x,i2,4x),/)                           
      aa=0.d0                                                           
      n1=nclass-1                                                       
      write(6,1540) aa,(aval(ii),ii=1,n1)                               
1540  format(12x,'inf:',10(2x,e8.2))                                    
      write(6,1550) (aval(ii),ii=1,nclass)                              
1550  format(12x,'sup:',10(2x,e8.2))                                    
      write(6,1030)                                                     
      write(6,1030)                                                     
      write(6,1560) (nom(ii),ii=1,nclass)                               
1560  format(1x,'ntotal',9x,10i10)                                      
      write(6,1030)                                                     
      write(6,1030)                                                     
      do 1600 kk=1,metat                                                
      write(6,1570) numero(kk),(nam(ii,kk),ii=1,nclass)                 
1570  format(1x,'etat ',i2,3x,'n:   ',10i10)                            
      if(.not.yen) then                                                 
      write(6,1580) (tabmp(ii,kk),ii=1,nclass)                          
      else                                                              
      write(6,1580) (tabmp(ii,kk),ii=1,nclass)                          
      write(6,1590) (tabvp(ii,kk),ii=1,nclass)                          
      endif                                                             
1580  format(11x,'emp: ',10(1x,f9.6))                                   
1590  format(11x,'evp: ',10(1x,f9.6))                                   
1600  write(6,1030)                                                     
      write(6,1030)                                                     
1800  continue                                                          
      write(6,1030)                                                     
      write(6,1700)                                                     
1700  format(1x,130('*'))                                               
      if(yprt) then                                                
      do 71 k=1,ncf                                                     
      l1=nd(k)+1                                                        
       l2=nd(k)+ne(k)                                                   
71     write(6,72) ne(k),nd(k),(trou(l),part(l),l=l1,l2)                
72    format(1x,30i4)                                                   
      do 73 i=1,mnorb                                                   
73    write(6,*) iorb(i),ispin(i),itsym(i),yocs(i)                      
      do 75 i=1,ntrsy                                                   
75    write(6,*) (yef(j,i),j=1,mnorb)                                   
      do 76 i=1,ncf                                                     
 76    write(6,*) (yocd(j,i),j=1,mnorb)                                 
      end if
      
      call fdate(timeday)
      write(6,9998) timeday
 9998 format(/,80(1h*),/,8(10h fin cipsi),/, 10x,a25,/,
     * 80(1h*))
      stop                                                              
      end                                                               
