      program fock
      implicit real*8 (a-h,o-x,z),logical*1 (y)                         
      character *26 timeday
      character*40 typener
      include 'pshf.prm'
      parameter (thr=1.d-7)
      logical chstan
      dimension f(doa,doa),ajk(doa,doa)                         
      integer*2 indic,jndic,nad,kt                                      
c    &,lndic                                                            
      common/nodim/norb,noca,nocb,igel,ncf,isym,nsym,isz,ntrsy          
      common/sy/ itsym(doa),its(nsymz,nsymz),mdegen(5),isydeg(nsymz,ns
     #ymz),itsyp(doa)                                                 
      common/truc/indic(doa*(doa+1)/2),jndic(doa*(doa+1)/2),
     &kt(10000),nbo(99),ndeb(500),nad(10000),
     &lndic(doa*(doa+1)/2)
      real*4 bijkl(kget)
      num(i)=i*(i-1)/2
      call openf

      call fdate(timeday)
      write(6,9999)timeday
      rewind 40                                                         
c            
      read (25) nsym,norb,noc,ntrans,(itsym(i),i=1,norb),noa,ydp,chstan
      write (40) nsym,norb,noc,ntrans,(itsym(i),i=1,norb),noa,ydp,chstan
      noca=noc
      nocb=noc
      ntyp=nsym                                                         
      ns=nsym   
       read(25) ((its(i,j),i=1,ns),j=1,ns)                              
       write(40) ((its(i,j),i=1,ns),j=1,ns)                             
      ntrsy=0                                                           
      if(ntrans.eq.0) go to 16                                          
      do 15 i=1,ntrans                                                  
      read (25) j,((isydeg(k,l+ntrsy),k=1,ntyp),l=1,j)                  
c      write(6,*) j,((isydeg(k,l+ntrsy),k=1,ntyp),l=1,j)                
      write(40) j,((isydeg(k,l+ntrsy),k=1,ntyp),l=1,j)                  
15    ntrsy=ntrsy+j                                                     
16    continue      
      nijkl=0
      if(chstan) then
        call reijkl_stan(nijkl) 
      else
	call reijkl_new(nijkl)
      end if
      iwant=nijkl                                                       
      write(6,910) iwant
910   format('0 taille maximum pour le bloc pqrs ou pqkl (variable)',i8)
      if(ydp) then
           if(kget.lt.iwant*2) then
           write(6,*) 'la dimension kget ',kget,' de bijkl doit etre 
     *     superieure a deux fois iwant( cas double precision) ',2*iwant
           stop
           end if
       else
           if(kget.lt.iwant) then
           write(6,*) 'la dimension kget ',kget,' de bijkl doit etre 
     *     superieure a iwant( cas simple precision) ',iwant
           stop
           end if
       endif
      if (chstan)then
         call itijkl_stan(bijkl,bijkl,ydp)
      else
	 call itijkl_new(bijkl,bijkl,ydp)
      end if
      write(*,*) 'test'
      read(25) enuc,ljkf    
      write(*,*) 'test',norb      
      read(50) ((f(i,j),j=1,i),i=1,norb)   
      write(*,*) 'test'      
      write(40) ((f(i,j),j=1,i),i=1,norb)                               
      if(ljkf.ne.0)then
      write(6,*)' integrales monoelectroniques'
      do i=1,norb
	write(6,711)i,(f(i,j),j=1,i)
      end do
      end if
      val=0.                                                            
      do 500 i=1,noc    
      val=val+f(i,i)*2.d0                                               
      tt=ai(i,i,i,i,bijkl,bijkl,ydp)         
      f(i,i)=f(i,i)+tt      
      val=val+tt                                                        
      im=i-1                                                            
      if(im.eq.0) go to 500                                             
      do 460 k=1,im              
      tt=ai(i,i,k,k,bijkl,bijkl,ydp)         
      tk=ai(i,k,i,k,bijkl,bijkl,ydp)         
      tt=tt+tt-tk                                                       
      f(i,i)=f(i,i)+tt                                                  
      f(k,k)=f(k,k)+tt 
      val=val+tt+tt  
      do 450 j=1,noc                                                    
450   f(i,k)=f(i,k)+ai(i,k,j,j,bijkl,bijkl,ydp)*2.d0-
     *ai(i,j,k,j,bijkl,bijkl,ydp)
460   continue
500   continue                                                          
      etot=val+enuc                                                     
      write(6,730) noc,val,etot                                         
      i1=noc+1                                                          
      do 600 i=i1,norb                                                  
      do 600 j=1,i                                                      
      tt=0.                                                             
      do 550 k=1,noc                                                    
550   tt=tt+2.d0*ai(i,j,k,k,bijkl,bijkl,ydp)-
     * ai(i,k,j,k,bijkl,bijkl,ydp)
600   f(i,j)=f(i,j)+tt                                                  
      do 605 i=1,norb
      do 605 j=1,i
      ajk(i,j)=ai(i,i,j,j,bijkl,bijkl,ydp)
      ajk(j,i)=ai(i,j,i,j,bijkl,bijkl,ydp)
605   continue
      if(ljkf.eq.0) go to 610                                           
      write(6,740)                                                      
740   format(' integrales j ' )                                         
      do 700 i=1,norb                                                   
  700 write(6,711) i,(ajk(i,j),j=1,i)                                   
      write(6,750)                                                      
  750 format(/,' integrales k')                                         
      do 705 i=1,norb                                                   
  705 write(6,711) i,(ajk(j,i),j=1,i)                                   
  610 continue                                                          
      do 3420 i=1,norb
      do 3420 j=1,i
      if(dabs(f(i,j)).lt.thr) f(i,j)=0.d0
 3420 continue
      write(6,760)                                                      
  760 format(/,' operateur de fock')                                    
      do 710 i=1,norb                                                   
  710 write(6,711) i,(f(i,j),j=1,i)                                     
  711 format(i5,12f10.6,(/,5x,12f10.6))                                 
730   format(//,130(1h*),/,' operateur de fock et energie totale',/,    
     &' ******attention******attention******attention******',/,         
     &' ces deux quantites sont calculees par rapport a un determinant c
     &onstruit sur des orbitales doublement occupees',/,                
     &      '   elles ne coincident avec les valeurs obtenues a la fin d
     &u calcul scf que dans le cas de couches fermees',/,               
     &' nombre de couches doublement occupees',i5,/,                    
     &' energie electronique',f20.8,/,                                  
     &' energie totale      ',f20.8,/,130(1h*))                         
      write(40) ((ajk(i,j),j=1,i),i=1,norb)                             
      write(40) ((ajk(j,i),j=1,i),i=1,norb)                             
c      transfert de f dans ajk pour probleme de precision               
      write(40) ((f(i,j),j=1,i),i=1,norb)                               
c   ecriture des integrales ijkl sur 40                                 
      if(chstan)then
	write (40) etot
        call wijkl_stan(bijkl,bijkl,ydp)
      else
	call wijkl_new(bijkl,bijkl,ydp)
      end if
      ij=num(norb)+norb                                                 
      typener='energie du determinant de reference'
      call stkener(typener,etot,1)
      call fdate(timeday)
       write(6,9998)timeday
 9999 format(/,80(1h*),/,8(10h   fock   ),/, 10x,a25,/,
     * 80(1h*))
 9998 format(/,80(1h*),/,8(10h fin fock ),/, 10x,a25,/,
     * 80(1h*))
      stop                                                              
      end                                                               
