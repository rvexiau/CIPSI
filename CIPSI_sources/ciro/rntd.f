      subroutine rntd(ii,jj,nbdif,mu,mv,isign)                          
      implicit real*8 (a-h,o-x,z),logical*1(y)                          
      include 'pshf.prm'
      common/ic/e(metz2),trou(8*ndimh),part(8*ndimh),
     & ne(2*ndimh),nd(2*ndimh),ispin(2*doa), 
     & iorb(2*doa),itsym(2*doa),its(nsymz,nsymz),igels(2),icount,
     & norb,noca,nocb,mnorb,ncf,nsym,isym,ntrsy,metat1,metat2,
     & yoc(2*doa),yprt              
      integer*4 ne,nd,iorb,ispin,itsym,its                    
      integer*4 trou,part                    
c        but                                                            
c                                                                       
c                  ssp de recherche des elements d interaction          
c                  entre deux determinants,numerote ii  et  jj          
c         pour l element de matrice <ii/r/jj>                           
c                                                                       
      dimension ita(20),itb(20),nv1(2),nv2(2),natu1(2),natu2(2)         
      equivalence(nv11,nv1(1)),(nv12,nv1(2)),(nv21,nv2(1)),(nv22,nv2(2))
c                                                                       
      i=ii                                                              
      j=jj                                                              
c                                                                       
c                                                                       
c                  le determinant - 1 - est choisi comme le             
c                  determinant ayant le plus d o.m. excitees            
c                                                                       
      if (ne(i)-ne(j) ) 5,5,7                                           
 7    kr = i                                                            
      kl= j                                                             
      yinv=.false.                                                      
      go to 9                                                           
 5    kr= j                                                             
      kl= i                                                             
      yinv=.true.                                                       
 9    ne1=ne(kr)                                                        
c                  test sur le fondamental                              
      if(ne1.eq.0) goto 44                                              
      ne2=ne(kl)                                                        
c                  nbdif : est le nombre de spin-orbitales              
c                          differentes                                  
      nbdif=ne1-ne2                                                     
c                                                                       
      if(nbdif.gt.1) goto 24                                            
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
      if(nbdif.gt.1) go to 24                                           
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
c      write(6,*)' ni,ne1,k-1=',ni,ne1,k-1                              
c      write(6,*)' ita(ni+ne1*(k-1))=',ita(ni+ne1*(k-1))                
c      write(6,*)' nv2(',i,')=',nv2(i)                                  
                                                                        
c                                                                       
      if(k.eq.2) go to 38                                               
40    nv1(i)=nv2(i)                                                     
      nv2(i)=nu                                                         
      ysig=.not.ysig                                                    
   38 continue                                                          
34    continue                                                          
c      write(6,*)' nv11,nv21=',nv11,nv21                                
      n1=iorb(nv21)                                                     
      n2=iorb(nv11)                                                     
      ns1=ispin(nv11)                                                   
      ns2=ispin(nv21)                                                   
      ncr=1                                                             
      if(ysig) ncr=-1                                                   
      isign=ncr                                                         
      if(yinv) go to 505                                                
      mu=n1                                                             
      mv=n2                                                             
      go to 24                                                          
505   mv=n1                                                             
      mu=n2                                                             
24    return                                                            
44    write(6,1000)                                                     
1000  format(/,1x,'erreur dans rntd:identite des determinants',/)       
      return                                                            
      end                                                               
