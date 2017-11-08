      subroutine ic001    
c         Replace subroutine given by diagonaliser for Lapack call; dynamic allocation for hmat 
      implicit real*8(a-h,o-x,z),logical*1 (y)      
      include 'pshf.prm'
      common/giv/nblock(ndetz),nstart(ndetz)                            
      common/spo/fdiag(2*doa),hdiag(ndetz),tocom(2*doa),coeff,      
     1ndef,ispin(2*doa),iwspin(2*doa),numero(metz),          
     2ne(ndetz+1),nd(ndetz+1),trou(ntpz),part(ntpz),iorb(2*doa),      
     3itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),yocs(2*doa),       
     4ygen(ndetz),yocd(2*doa,ndetz),yspin,yef(ndetz,nsymz)                   
      integer*4 ne,nd,trou,part,iorb,iwspin,itsym,its,isytr,nesp,nespo, 
     1kq,ner,ndr,tr,pr                                                  
      character*4 ibla2,ing,ig,mig                                                                
      common/cipsi/escf,test,teseff,tdet,twdet,ndiag,norb,mnorb,
     1             noca,nocb,noc1,noc2,igela,igelb,ncf,nsym,isym,isz,
     2             ntrsy,metat,lecdet,maxdet,yprt,yion,ypertu,ybrd,yen,
     3             ywvi,ywvf,ytoul,ystkic,ywdeg
      common/ic/c(metz*ndetz),e(ndetz),eb(metz),e2b(metz),emp(metz),    
     1teste(metz),e2mp(metz),dpivot,ietat,ishif,e2vp(metz),             
     2fmpb(2*doa),hefmp(nhefz),hefvp(nhefz),hefb(nhefz),              
     3psimp(metz),psivp(metz)                                           
      common/infogv/nbdet(9),iord60(ndetz),ynoic                        
      dimension fmul(metz)                                              
c      dimension hic(ndetz),cx(nblokz),hmat(ndetz*(ndetz+1)/2)
      dimension hic(ndetz),cx(nblokz)      
      real*8,dimension(:),allocatable :: hmat
      data ibla1,ibla2/0,'=>'/                                          
      data ing,ig/'ng',' g'/      
      write(6,5888)                                                     
      write(6,1300)                                                     
      write(6,1310)                                                     
1310  format(1x,'ic variationnelle')                                    
      write(6,1300)                                                     
1300  format(1x,17('*'))                                                
      write(6,5888)                                                     
5888  format(/)                                                         
      if(ncf.gt.ndetz) then                                             
      write(6,5556) ndetz,ncf                                           
5556  format(1x,'dimension maximale du sous espace variationnel:ndetz=' 
     *,i9,' ncf=',i9)                                                   
      stop                                                              
      endif              
      if (ynoic) then                                                   
c--- lecture sur f59 des vp et valp issues de bdav, puis reclassement   
      rewind 4                                                          
      read(4) n,n,n,netat,mcf                                           
         rewind 4    

         rewind 59                                                      
         read(59) ncf,metat,(c(i),i=1,metat*ncf),(e(i),i=1,metat)       
         write(6,*) ' resultats variationnels lus sur file 59'    
         write(6,*) ' nombre d etats ',metat   ,ncf       
         e(1:metat)=e(1:metat)-escf
         letat=metat  
       if(ncf.ne.mcf) then                                              
         write(6,*)                                                     
            write(6,*) ' files 4 et 59 incompatibles'                   
            write(6,*) ' file 04 ncf: ',mcf                             
            write(6,*) ' file 59 ncf: ',ncf                             
            stop  
         endif        
       letat=metat                                                     
       do 5594 i=1,ncf                                                  
5594   hdiag(i)=hdig(i)                                                 
       else                                                             
                                                                        
c                                                                       
c   construction et stockage sur file 1 de la matrice adiabatique       
c                  
      allocate(hmat(ncf*(ncf+1)/2))
      imax=0                                                            
      nsize=nblokz                                                      
      ng=1                                                              
      ij=0
      do 350 i=1,ncf             
      hic(i)=hdig(i)                                                    
      hdiag(i)=hic(i)                                                   
      ij=ij+1
      hmat(ij)=hdiag(i)
      if(i.eq.ncf) go to 310                                            
      ip=i+1           
      do 304 j=ip,ncf         
      hic(j)=0.d0                                                       
      ndj=nd(j)+1                                                       
      nej=ne(j)                                                         
      ndif=ne(i)-nej                                                    
      nej=nej+ndj-1                                                     
      do 300 k=ndj,nej                                                  
      ntp=trou(k)                                                       
      if(.not.yocd(ntp,i))ndif=ndif+1                                   
      ntp=part(k)                                                       
      if(.not.yocd(ntp,i))ndif=ndif+1                                   
300   continue          
      if(ndif.le.2) hic(j)=hntd(i,j)     
      ij=ij+1 
      hmat(ij)=hic(j)
304   continue                                                          
310   continue
      if(ystkic) write(6,3222) (hic(j),j=i,ncf)                         
3222  format(10(1x,f10.6))                                              
350   continue                                                          
c                                                                       
c    si ystkic=t, sauvegarde de la matrice d'ic sur 1                   
c    sans diagonalisation                                               
      if(.not.ystkic) go to 3836                                        
      write(6,3839) (hdiag(i),i=1,ncf)                                  
 3839 format(/,'  elements diagonaux ',/,(10f12.6))                     
      stop                                                              
 3836 continue                                                          
c   diagonalisation                                                     
c  *****************                                                    
c                                                                       
580   if(metat.gt.ncf) metat=ncf                                        
      nme=numero(metat)                                                 
      maxe=max0(nme,ndiag)                                              
      if(maxe.gt.ncf) maxe=ncf                                          
      maxe=maxe-1                                                       
      letat=min0(ncf,((maxe/6)*6+6),((ndetz*metz)/ncf))   
!      call given(hmat,e,c,ncf,ncf,letat) 
      call diagonaliser(ncf,hmat,e(1:ncf),c(1:ncf*letat),letat) 
      deallocate(hmat)
      endif            
      write(6,5777) ncf                                                 
5777  format(1x,'nombre de determinants   ncf=',i9)                     
      write(6,5888)                                                     
      write(6,*) ' valeurs propres'                                     
      write(6,103)(e(i),i=1,3*metat)                                        
      if(.not.ywvi) go to 111                                           
c                                                                       
c   impression des vecteurs initiaux (ywvi=f,defaut)                    
c                       
      max=0                                                             
29    max1=max+1                                                        
      max=max1+5                                                        
      if(letat.lt.max) max=letat                                        
      write(6,1006)                                                     
      do 12 k=1,ncf                                                     
      mono1=1+nd(k)                                                     
      mono2=mono1+ne(k)-1                                               
      write(6,100) (c(ncf*(j-1)+k),j=max1,max)                          
      write(6,6666) (iorb(trou(i)),iwspin(trou(i)),i=mono1,mono2),      
     * ibla1,ibla2, (iorb(part(i)),iwspin(part(i)),i=mono1,mono2)       
12    continue                                                          
      if(max.lt.letat) go to 29                                         
111   continue                                                          
c                                                                       
c   symmetrisation des vecteurs                                         
c                                                                       
      call degen(letat)                   

      nbr=metat*(metat+1)/2                                             
      write(6,5888)                                                     
      write(6,5888)                                                     
      write(6,1420)                                                     
      write(6,1410)                                                     
1410  format(1x,'caracteristiques d ordre zero des etats ',             
     *'traites en perturbation')                                        
      write(6,1420)                                                     
1420  format(1x,63('*'))                                                
      write(6,*)                                                        
      write(6,3776) twdet                                               
3776  format(1x,'coefficients imprimes > twdet=',f5.3,/)                
      if(.not.ywvf) go to 740                                           
c                                                                       
c    impression des vecteurs perturbes apres symmetrisation             
c    et/ou reclassement ( ywvf=t,defaut)                                
c                                                                       
      write(6,5888)                                                     
      write(6,1100)                                                     
1100  format(1x,//,' energies: ')                                       
      write(6,103) (e(i),i=1,metat)                                     
      max=0                                                             
129   max1=max+1                                                        
      max=max1+5                                                        
      if(metat.lt.max) max=metat                                        
      write(6,1200)                                                     
1200  format(1x,'vecteurs propres utiles (et symmetrises)',/)           
      do 112 k=1,ncf 
      ygen(k)=.true.
      mncf=-ncf                                                         
      do 6820 ii=1,metat                                                
      mncf=mncf+ncf                                                     
      if(dabs(c(mncf+k)).gt.twdet) go to 6822                           
6820  continue                                                          
      go to 112                                                         
6822  continue                                                                                                         
      mono1=1+nd(k)                                                     
      mono2=mono1+ne(k)-1                                               
      write(6,1400) ig,(c(ncf*(j-1)+k),j=max1,max)                     
      write(6,6666) (iorb(trou(i)),iwspin(trou(i)),i=mono1,mono2),      
     * ibla1,ibla2, (iorb(part(i)),iwspin(part(i)),i=mono1,mono2)       
112   continue                                                          
      if(max.lt.metat) go to 129                                        
      write(6,5888)                                                     
      write(6,5888)                                                     
6666  format(1x, 8(i3,a2,' '),/,(1x,8(i3,a2,' ')))                      
100   format(49x,3x,6f12.8)                                             
103   format(56x,6f12.8)                                                
1400  format(50x,a2,4x,6f12.8)                                          
101     format('0',(47x,6i14))                                          
1005  format(1x,'energies variationnelles :')                           
1006  format(/,1x,'vecteurs: '/)                                        
740   continue                                                          
      return                                                            
      end                                                               
