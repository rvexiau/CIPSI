      function hntd(ii,jj)                                              
      implicit real*8 (a-h,o-x,z),logical*1 (y)                         
      include 'pshf.prm'
      common/jkf/ f(doa*(doa+1)/2),aj(doa*(doa+1)/2),           
     2ak(doa*(doa+1)/2)                                             
c        but                                                            
c                                                                       
c                  ssp de recherche des elements d interaction          
c                  entre deux determinants,numerote ii  et  jj          
c                                                                       
      dimension ita(20),itb(20),nv1(2),nv2(2),natu1(2),natu2(2)         
      equivalence(nv11,nv1(1)),(nv12,nv1(2)),(nv21,nv2(1)),(nv22,nv2(2))
      common/spo/fdiag(2*doa),hdiag(ndetz),tocom(2*doa),coeff,      
     1ndef,ispin(2*doa),iwspin(2*doa),numero(metz),          
     2ne(ndetz+1),nd(ndetz+1),trou(ntpz),part(ntpz),iorb(2*doa),      
     3itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),yocs(2*doa),       
     4ygen(ndetz),yocd(2*doa,ndetz),yspin,yef(ndetz,nsymz)                   
      common/degene/nesp(ndetz),nespo(ndetz),isytr(ndetz,nsymz),kkbrd,  
     1ysaut(ndetz)                                                      
      integer*4 ne,nd,trou,part,iorb,iwspin,itsym,its,isytr,nesp,nespo  
      common/cipsi/escf,test,teseff,tdet,twdet,ndiag,norb,mnorb,
     1             noca,nocb,noc1,noc2,igela,igelb,ncf,nsym,isym,isz,
     2             ntrsy,metat,lecdet,maxdet,yprt,yion,ypertu,ybrd,yen,
     3             ywvi,ywvf,ytoul,ystkic,ywdeg
      common/ic/c(metz*ndetz),e(ndetz),eb(metz),e2b(metz),emp(metz),    
     1teste(metz),e2mp(metz),dpivot,ietat,ishif,e2vp(metz),             
     2fmpb(2*doa),hefmp(nhefz),hefvp(nhefz),hefb(nhefz),              
     3psimp(metz),psivp(metz)                                           
      common/truc/indic(doa*(doa+1)/2),jndic(doa*(doa+1)/2),    
     1ndeb(500),nad(2000),kt(2000),                                     
     2nbo(99)                                                           
      integer*4 indic,jndic,nad,kt                                      
      logical iden,iprin1      
      num(i)=i*(i-1)/2                
      i=ii                                                              
      j=jj                                                              
      iden=.false.                                                      
c                                                                       
c                  le determinant - 1 - est choisi comme le             
c                  determinant ayant le plus d o.m. excitees            
c                           
      if (ne(i)-ne(j) ) 5,5,7                                           
 7    kr = i                                                            
      kl= j                                                             
      go to 9                                                           
 5    kr= j                                                             
      kl= i                                                             
 9    ne1=ne(kr)                                                        
c                  test sur le fondamental                              
      if(ne1.eq.0) goto 44                                              
      ne2=ne(kl)                                                        
c                  nbdif : est le nombre de spin-orbitales              
c                          differentes                                  
      nbdif=ne1-ne2                                                     
c                                                                       
      if(nbdif.gt.2) goto 24                                            
c                  construction pour le determinant -1- et -2-          
c                  de it  et ns,qui contiennent les numeros             
c                  des orbitales et leurs spins                         
c                                                                       
      nd1=nd(kr)                                                        
      nd2=nd(kl)                                                        
      do 15 j=1,ne1                                                     
      ita(j)=trou(j+nd1)                                                
      ita(j+ne1)=part(j+nd1)                                            
15    continue                            
c                                                                       
c                  ncr :  est le nombre de croisements                  
c                         dans le diagramme d interaction               
      ysig=.false.                                                      
c                  si ne2=0 le determinant -2- est le fondamental       
      if (ne2.eq.0) goto 125                                            
121   continue                                                          
      do 25 j= 1,ne2                                                    
      itb(j)=trou(j+nd2)                                                
25    itb(j+ne2)=part(j+nd2)                                            
    4 do 8 k=1,2                                                        
      k1=(k-1)*ne1                                                      
      k2=(k-1)*ne2                                                      
      do 10i=1,ne2                                                      
c                                                                       
      nit1=itb(k2+i)                                                    
      j1=1+k1                                                           
      j2=ne1+k1                                                         
      do 12 j=j1,j2                                                     
      njt2=ita(j)                                                       
      if(nit1.ne.njt2) go to 12                                         
      if(i.eq.(j-k1)) go to 10                                          
      ita(j)=ita(i+k1)                                                  
      ita(i+k1)=njt2                                                    
      ysig=.not.ysig                                                    
      go to 10                                                          
   12 continue                                                          
      nbdif=nbdif+1                                                     
      if(nbdif.gt.2) go to 24                                           
      nv1(nbdif)=nit1                                                   
c                                                                       
      natu1(nbdif)=k                                                    
      natu2(nbdif)=i                                                    
   10 continue                                                          
    8 continue                                                          
125   if(nbdif.le.0) go to 44                                           
c                                                                       
c                   nar : excitation suplementaire du determinant -1-   
c                        par rapport au determinant -2-                 
   26 nar=ne1-ne2                                                       
      if(nar.le.0) go to 28                                             
   30 do 32 i=1,nar                                                     
      net=ne1+1-i                                                       
      nv1(i)=ita(net)                                                   
      nv2(i)=ita(net+ne1)                                               
   32 continue                                                          
   28 nbar=nbdif-nar                                                    
      if(nbar.le.0) go to 34                                            
   36 nar1=nar+1                                                        
      do 38 i=nar1,nbdif                                                
      k=natu1(i)                                                        
      ni=natu2(i)                                                       
      nu=nv1(i)                                                         
      nv2(i)=ita(ni+ne1*(k-1))                                          
c                                                                       
      if(k.eq.2) go to 38                                               
40    nv1(i)=nv2(i)                                                     
      nv2(i)=nu                                                         
      ysig=.not.ysig                                                    
   38 continue                                                          
34    continue                                                          
c                  calcul de l element de matrice                       
 2    ha=0.                                                             
      n1=iorb(nv11)                                                     
      n2=iorb(nv21)                                                     
      ns1=ispin(nv11)                                                   
      ns2=ispin(nv21)                                                   
      ncr=1                                                             
      if(ysig) ncr=-1                                                   
      if(nbdif.eq.1) go to 58                                           
c                  les -2- determinants differe par -2- spin-orbitales  
60    n3=iorb(nv12)                                                     
      n4=iorb(nv22)                                                     
      ns3=ispin(nv12)                                                   
      ns4=ispin(nv22)           
      if(ns1.eq.ns2)ha=ai(n1,n2,n3,n4)                                  
      if(ns1.eq.ns4)ha=ha-ai(n1,n4,n3,n2) 
      go to 75                                                          
c                                                                       
c                  les -2- determinants differe par -1- spin-orbitale   
58    if(n2.ge.n1) go to 57                                             
      n12=num(n1)+n2                                                    
      go to 56                                                          
57    n12=num(n2)+n1                                                    
56    ha=f(n12)                                                         
      do 72 i=1,ne1                                                     
      np=iorb(ita(i))                                                   
      nsu=ispin(ita(i))                                                 
      aix=ai(n1,n2,np,np)                                             
      ha=ha-aix
      if(ns1.eq.nsu) ha=ha+ai(n1,np,n2,np)                              
      np=iorb(ita(i+ne1))                                               
      if(np.gt.norb) go to 72                                           
      nsu=ispin(ita(i+ne1))                                             
      ha=ha+ai(n1,n2,np,np)                                             
      if(ns1.eq.nsu) ha=ha-ai(n1,np,n2,np)                              
   72 continue                                                          
75    hntd=ha*ncr                          
      if(.not.yprt) return                                              
      write(6,1007) hntd                                                
 1007 format(10x,f15.10)                                                
      return                                                            
 24   hntd=0.0                                                          
      if(.not.yprt) return                                              
      write(6,1005)                                                     
 1005 format(1h0,'pas d''interaction')                                  
      return                                                            
   44 iden=.true.                                                       
      hntd=0.0                                                          
      if(.not.yprt) return                                              
      write(6,1006)                                                     
 1006 format (1h0,'identite')                                           
      return                                                            
      end                                                               
