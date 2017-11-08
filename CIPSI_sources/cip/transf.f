      subroutine transf(nec)                                            
      implicit real*8(a-h,o-x,z),logical*1(y)                           
      include 'pshf.prm'
c   cette subroutine teste la validite des excitations                  
c  elle elimine les determinants non orthodoxes                         
c  et ne conserve que les excitations conformes                         
c                                                                       
      integer*4 ne,nd,trou,part,iorb,iwspin,itsym,its,isytr,nesp,nespo  
      common/cipsi/escf,test,teseff,tdet,twdet,ndiag,norb,mnorb,
     1             noca,nocb,noc1,noc2,igela,igelb,ncf,nsym,isym,isz,
     2             ntrsy,metat,lecdet,maxdet,yprt,yion,ypertu,ybrd,yen,
     3             ywvi,ywvf,ytoul,ystkic,ywdeg
      common/spo/fdiag(2*doa),hdiag(ndetz),tocom(2*doa),coeff,      
     1ndef,ispin(2*doa),iwspin(2*doa),numero(metz),          
     2ne(ndetz+1),nd(ndetz+1),trou(ntpz),part(ntpz),iorb(2*doa),      
     3itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),yocs(2*doa),       
     4ygen(ndetz),yocd(2*doa,ndetz),yspin                                    
      common/degene/nesp(ndetz),nespo(ndetz),isytr(ndetz,nsymz),kkbrd,  
     *ysaut(ndetz),imul(nmulz),yef(ndetz,nsymz)                         
      common     kf,km,kd,kt,kq,kp,kh,k7,k0,nt(9),np(9),mt(9),mp(9)     
      equivalence (k,ncf)                                               
      mnorb=norb+norb                                                   
      ncfi=k+1                                                          
      if(nec-1) 10,30,20                                                
10    k=k+1                                                             
      ne(k)=0                                                           
      nd(k)=0                                                           
      trou(1)=0                                                         
      part(1)=0                                                         
      if(isz.ne.0) k=k-1                                                
      isym=1                                                            
      return                                                            
c     les spinorbitales alpha sont numerotees de 1 a norb               
c     les spinorbitales beta de norb+1 a mnorb=2*norb                   
c     ordonnons trous et particulesbpar numeros croissants              
20    necm=nec-1                                                        
      do 25 i=1,necm                                                    
      ip=i+1                                                            
      do 25 j=ip,nec                                                    
      if(nt(j).gt.nt(i)) go to 24                                       
      n=nt(j)                                                           
      nt(j)=nt(i)                                                       
      nt(i)=n                                                           
24    if(np(j).gt.np(i)) go to 25                                       
      n=np(j)                                                           
      np(j)=np(i)                                                       
      np(i)=n                                                           
25    continue                                                          
30    is=1                                                              
      ls=0                                                              
c                                                                       
c  test sur les numeros des trous et des particules                     
c                                                                       
      if(nt(1).lt.igela) go to 80                                       
      if(np(nec).gt.mnorb.and..not.yion)go to 80                        
      if(nt(nec).gt.nocb) go to 80                                      
      if(np(nec).le. noca) go to 80                                     
      do 35 i=1,nec                                                     
      if(nt(i).le.noca) go to 31                                        
      if(nt(i).lt.igelb) go to 80                                       
31    if(np(i).gt.nocb) go to 32                                        
      if(np(i).gt.norb) go to 79                                        
      if(np(i).le.noca) go to 80                                        
c                                                                       
c test sur la symetrie et le spin des excitations                       
c                                                                       
32    is=its(is,itsym(nt(i)))                                           
      is=its(is,itsym(np(i)))                                           
      ls=ls+ispin(np(i))                                                
      ls=ls-ispin(nt(i))                                                
35    continue                                                          
      if(isym.eq.0) isym=is                                             
      if(.not.yion.and.is.ne.isym) go to 90                             
      if(ls.ne.isz) go to 96                                            
c  on garde les bons                                                    
      if(k.eq.0) go to 50                                               
      do 40 i=1,k                                                       
      if(ne(i).ne.nec) go to 40                                         
      ndk=nd(i)                                                         
      do 39 j=1,nec                                                     
      if(trou(ndk+j).ne.nt(j)) go to 40                                 
      if(part(ndk+j).ne.np(j)) go to 40                                 
39    continue                                                          
      go to 70                                                          
40    continue                                                          
      ndk=nd(k)+ne(k)                                                   
      if(ne(k).eq.0) ndk=nd(k)+1                                        
      go to 51                                                          
50    ndk=0                                                             
51    k=k+1                                                             
      if(k.gt.ndetz) go to 56                                           
      do 55 j=1,mnorb                                                   
55    yocd(j,k)=.false.                                                 
      ne(k)=nec                                                         
      nd(k)=ndk                                                         
      do 60 j=1,nec                                                     
c                                                                       
      yocs(nt(j))=.true.                                                
      yocs(np(j))=.true.                                                
      yocd(nt(j),k)=.true.                                              
      yocd(np(j),k)=.true.                                              
59    trou(ndk+j)=nt(j)                                                 
60    part(ndk+j)=np(j)                                                 
c                                                                       
c                                                                       
70    continue                                                          
120   return                                                            
79    if(yion) return                                                   
80    write(6,85)                                                       
c                                                                       
85    format(' 1 spin orbitale est hors de l espace de base')           
      write(6,*)(nt(i),np(i),i=1,nec)                                   
      stop                                                              
90    write(6,95)                                                       
      write(6,*) (nt(i),np(i),i=1,nec)                                  
95    format(' dterminant de mauvaise symetrie')                        
      stop                                                              
96    write(6,97)                                                       
97    format(' mauvais spin')                                           
      write(6,*)(nt(i),np(i),i=1,nec)                                   
      stop                                                              
56    write(6,99)  ndetz                                                
                                                                        
99    format(1x,'trop de determinants dans s: maximum ',i5,/)           
      stop                                                              
      end                                                               
