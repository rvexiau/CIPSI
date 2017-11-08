      subroutine spasym(netat)                                          
      implicit real*8(a-h,o-x,z),logical*1(y)                           
      include 'pshf.prm'
      common/degene/nesp(ndetz),nespo(ndetz),isytr(ndetz,nsymz),kkbrd,  
     1ysaut(ndetz),imul(nmulz),yheff(metz*(metz+1)/2),yef(ndetz,nsymz)  
      integer*4 isytr,nesp,nespo,imul,itm,itk,nti,ntl,npl               
      common/cipsi/escf,test,teseff,tdet,twdet,ndiag,norb,mnorb,
     1             noca,nocb,noc1,noc2,igela,igelb,ncf,nsym,isym,isz,
     2             ntrsy,metat,lecdet,maxdet,yprt,yion,ypertu,ybrd,yen,
     3             ywvi,ywvf,ytoul,ystkic,ywdeg
      common/ic/c(metz*ndetz),e(ndetz),eb(metz),e2b(metz),emp(metz),    
     1teste(metz),e2mp(metz),dpivot,ietat,ishif,e2vp(metz),             
     2fmpb(2*doa),hefmp(nhefz),hefvp(nhefz),hefb(nhefz),              
     3psimp(metz),psivp(metz)                                           
      common/spo/fdiag(2*doa),hdiag(ndetz),tocom(2*doa),coeff,      
     1ndef,ispin(2*doa),iwspin(2*doa),numero(metz),          
     2ne(ndetz+1),nd(ndetz+1),trou(ntpz),part(ntpz),iorb(2*doa),      
     3itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),yocs(2*doa),       
     4ygen(ndetz),yocd(2*doa,ndetz),yspin                                    
      common/nature/itm(metz,nsymz)                                     
      integer*4 ne,nd,trou,part,iorb,iwspin,itsym,its,npk,ntk           
      dimension nti(nsymz),itk(ndetz,10)                                
      dimension ntk(10),npk(10),ntl(10),npl(10)                         
      equivalence (itk(1,1),imul(1))                                    
      equivalence (itm(1,1),ntk(1)),(itm(11,1),npk(1))                  
      equivalence (ntl(1),nti(1)),(npl(1),nti(11))                      
c                                                                       
c   determination des transformes itk(k,l) des determinants k de s par  
c   les operations de symetrie l (numeros et signes)                    
c                                                                       
      isym=1                                                            
      nek=ne(1)                                                         
      ndk=nd(1)                                                         
        do 7200 i=1,nek                                                 
      istrou=itsym(trou(ndk+i))                                         
        isym=its(isym,istrou)                                           
      ipart=part(ndk+i)                                                 
      if(ipart.gt.mnorb) go to 7200                                     
      ispart=itsym(ipart)                                               
      isym=its(isym,ispart)                                             
7200  continue                                                          
      nutil=0                                                           
      if(ntrsy.le.1) return                                             
      if (nsym.eq.5) then                                               
      if (isym.eq.1) nutil=2                                            
      if (isym.eq.5) nutil=ntrsy+1                                      
      endif                                                             
      if (nsym.eq.10) then                                              
      if ((isym.eq.1).or.(isym.eq.2)) nutil=2                           
      if ((isym.eq.9).or.(isym.eq.10)) nutil=ntrsy+1                    
      endif                                                             
c      write(6,*) 'nsym,ntrsy,isym,mnorb',nsym,ntrsy,isym,mnorb         
      if(nutil.eq.ntrsy+1) then                                         
      do 8800 k=1,ncf                                                   
8800  yef(k,ntrsy+1)=.true.                                             
      do 8805 l=1,mnorb                                                 
         if (isytr(l,ntrsy+1).lt.0) then                                
            isytr(l,ntrsy+1)=-isytr(l,2)                                
         else                                                           
            isytr(l,ntrsy+1)=isytr(l,2)                                 
         endif                                                          
8805  continue                                                          
c      write(6,*) 'isytr'                                               
c      do 8810 ll=1,nutil                                               
c      write(6,*) (yef(k,ll),k=1,ncf)                                   
c8810  write(6,*) (isytr(l,ll),l=1,mnorb)                               
      endif                                                             
      yisz0=isz.eq.0                                                    
      if(.not.yisz0) l1=2                                               
      l1=ntrsy+1                                                        
      if(yion) isytr(mnorb+1,nutil)=mnorb+1                             
      do 100 k=1,ncf                                                    
      if(ysaut(k)) go to 100                                            
      nespk=nesp(k)+1                                                   
      nec=ne(k)                                                         
      if(nec.ne.0) go to 105                                            
      do 195 l=l1,nutil                                                 
195   itk(k,l)=k                                                        
      go to 100                                                         
105   nec1=nec-1                                                        
c      write(6,*) ' det ',k,' deg',nespk                                
      do 110 kg=1,nespk                                                 
      kk=k+kg-1                                                         
      ndk=nd(kk)                                                        
c      write(6,*) kk                                                    
      do 175 i=1,nec                                                    
      ntk(i)=trou(i+ndk)                                                
175   npk(i)=part(i+ndk)                                                
      do 120 l=l1,nutil                                                 
      isigne=1                                                          
      ncr=1                                                             
      do 125 i=1,nec                                                    
      nl=ntk(i)                                                         
      np=npk(i)                                                         
      ntli=isytr(nl,l)                                                  
      npli=isytr(np,l)                                                  
c      write(6,*) 'trous ',nl,ntli                                      
c      write(6,*) ' parts ',np,npli                                     
      if((ntli.ne.0).and.(npli.ne.0)) go to 155                         
      itk(kk,l)=0                                                       
      go to 120                                                         
155   continue                                                          
      if(ntli.le.0) isigne=-isigne                                      
      if(npli.le.0) isigne=-isigne                                      
      ntl(i)= abs(ntli)                                                 
      npl(i)= abs(npli)                                                 
125   continue                                                          
      if(yef(k,l)) go to 145                                            
       ii=kk                                                            
      go to 135                                                         
145   do 130 ig=1,nespk                                                 
      ii=k+ig-1                                                         
      do 140 i=1,nec                                                    
      it=ntl(i)                                                         
      ip=npl(i)                                                         
      if(.not.yocd(it,ii)) go to 130                                    
      if(.not.yocd(ip,ii)) go to 130                                    
140   continue                                                          
      go to 135                                                         
130   continue                                                          
      if(nec.eq.1) go to 198                                            
135   do 160 i=1,nec1                                                   
      ip=i+1                                                            
      do 150 j=ip,nec                                                   
      if(ntl(j).gt.ntl(i)) go to 165                                    
      n=ntl(j)                                                          
      ntl(j)=ntl(i)                                                     
      ntl(i)=n                                                          
      ncr=-ncr                                                          
165   if(npl(j).gt.npl(i)) go to 150                                    
      n=npl(j)                                                          
      npl(j)=npl(i)                                                     
      npl(i)=n                                                          
      ncr=-ncr                                                          
150   continue                                                          
160   continue                                                          
198   if(isigne.lt.0) ii=-ii                                            
      if(ncr.lt.0) ii=-ii                                               
      itk(kk,l)=ii                                                      
120   continue                                                          
110   continue                                                          
100   continue                                                          
      if(.not.yprt) go to 170                                           
      write(6,*) 'itk'                                                  
      write(6,180) (k,k=1,ncf)                                          
      do 185 l=l1,nutil                                                 
185   write(6,180) (itk(k,l),k=1,ncf)                                   
180   format(20(1x,i4))                                                 
170   continue                                                          
c                                                                       
c  determination du signe itm(m,l) des transformes des etats m          
c  par les operations de symetrie l                                     
c                                                                       
      do 200 m=1,netat                                                  
      mm=m                                                              
      mncf=(mm-1)*ncf                                                   
      do 210 l=l1,nutil                                                 
      x=0.d0                                                            
      do 220 k=1,ncf                                                    
      if(itk(k,l).eq.0) go to 220                                       
      kk=itk(k,l)                                                       
      ikk= abs(kk)                                                      
      ckk=c(mncf+k)*c(mncf+ikk)                                         
      if(kk.lt.0) ckk=-ckk                                              
      x=x+ckk                                                           
220   continue                                                          
      itm(m,l)=1                                                        
      if(x.lt.0.d0) itm(m,l)=-1                                         
210   continue                                                          
200   continue                                                          
c     if(.not.yprt) go to 250                                           
      do 230 l=l1,nutil                                                 
230   write(6,180) (itm(m,l),m=1,netat)                                 
250   continue                                                          
c                                                                       
c      if(nutil.eq.0) then                                              
c         write(6,*) ' etats pi'                                        
c      else                                                             
c         write(6,*) 'etats sigma ou delta -- itm'                      
c         write(6,180) (itm(m,nutil),m=1,netat)                         
c      endif                                                            
      return                                                            
      end                                                               
