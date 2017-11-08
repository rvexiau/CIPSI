      subroutine reijkl_stan(nijkl)                                     
      implicit real*8(a-h,o-x,z),logical*1(y)                           
      include 'pshf.prm'                                  
      common/truc/indic(doa*(doa+1)/2),jndic(doa*(doa+1)/2),
     & kt(10000),nbo(99),ndeb(500),nad(10000),
     & lndic(doa*(doa+1)/2)                                   
      integer*2 indic,jndic,nad,kt
      common/nodim/ norb,noca,nocb,igel,ncf,isym,nsym,isz,ntrsy         
      common nos(nsymz),ns(doa),nus(nsymz),nbp(99),npair,ipair(2000)   
     &,noms(doa,nsymz)                                                
      integer*2 noms                                                    
      common/sy/ itsym(doa),its(nsymz,nsymz),mdegen(5),               
     #isydeg(nsymz,nsymz),itsyp(doa)                                  
      num(i)=i*(i-1)/2
      read(25) npair,(ipair(i+i-1),ipair(i+i),i=1,npair)                
      write(40) npair,(ipair(i+i-1),ipair(i+i),i=1,npair)               
      read (25) nsyt,(nad(i),i=1,nsyt)                                  
      write(40) nsyt,(nad(i),i=1,nsyt)                                  
      read(25) (kt(i),i=1,nsyt)                                         
      write(40) (kt(i),i=1,nsyt)                                        
      ij=0                                                              
      kls=0                                                             
      do 9 lsym=1,nsym                                                  
      k=0                                                               
      m=0                                                               
      do 4 i=1,norb                                                     
      if(itsym(i).ne.lsym) go to 4                                      
      k=k+1                                                             
      noms(k,lsym)=i                                                    
c     if(.not.yocs(i).or..not.yocs(i+norb)) go to 4                     
      m=m+1                                                             
      ns(i)=m                                                           
4     continue                                                          
      nos(lsym)=k                                                       
      nus(lsym)=m                                                       
c k nombre total m nombre a occupation variable                         
      do 9 ksym=1,lsym                                                  
      kls=kls+1                                                         
      if(lsym.eq.ksym) go to 5                                          
      nbp(kls)=k*nos(ksym)                                              
      nbo(kls)=m*nus(ksym)+(k-m)*nus(ksym)+(nos(ksym)-nus(ksym))*m      
      go to 9                                                           
5     nbp(kls)=num(k)+k                                                 
      nbo(kls)=num(m)+m+m*(k-m)                                         
9     continue                                                          
72    format('  nos,  nus   ns   nbo   nbp   ipair')                    
c     write(2,*) 'nos:'
c     write(2,70) nos                                                   
c     write(2,*) 'nus:'
c     write(2,70) nus                                                   
c     write(2,*) 'ns:'
c     write(2,70) ns                                                    
c     write(2,*) 'nbo:'
c     write(2,70) nbo                                                   
c     write(2,*) 'nbp:'
c     write(2,70) nbp                                                   
c     write(2,*) 'ipair(npair=',npair,')'
c     write(2,70) (ipair(i),i=1,npair),(ipair(npair+i),i=1,npair)       
      nijkl=0                                                           
c     le premier indice dans la paire est l indice exterieur (colonne)  
      do 30 n=1,npair                                                   
      i=ipair(n+n-1)                                                    
      ni=nbo(i)                                                         
      j=ipair(n+n)                                                      
      nj=nbo(j)                                                         
      nij=ni*nj                                                         
      if(i-j) 15,17,19                                                  
15    ijkls=j*(j-1)/2+i                                                 
      go to 20                                                          
17    ijkls=i*(i+1)/2                                                   
      nadi=nad(ijkls)                                                   
      if(nadi.lt.0) nadi=-nadi                                          
c                                                                       
c                                                                       
      go to 21                                                          
19    ijkls=i*(i-1)/2+j                                                 
   20 continue                                                          
      nadi=nad(ijkls)                                                   
      if(nadi.lt.0) nadi=-nadi                                          
   21 if(kt(ijkls).eq.1.or.kt(ijkls).eq.10) nij=(nij+ni)/2              
      ndeb(nadi)=nijkl                                                  
      nijkl=nijkl+nij                                                   
30    continue                                                          
      ijs=0                                                             
      do 51 is=1,nsym                                                   
      ks=nos(is)                                                        
      do 51 js=1,is                                                     
      ijs=ijs+1                                                         
      ls=nos(js)                                                        
      if(ks.eq.0.or.ls.eq.0) go to 51                                   
      ijn=0                                                             
      do 50 k=1,ks                                                      
      i=noms(k,is)                                                      
      if(is.eq.js) ls=k                                                 
      do 50 l=1,ls                                                      
      j=noms(l,js)                                                      
      ijn=ls*(k-1)+l                                                    
      if(is.eq.js) ijn=num(k)+l                                         
      ijm=ks*(l-1)+k                                                    
      if(is.eq.js) ijm=ijn                                              
      ij=num(i)+j                                                       
      if(j.gt.i) ij=num(j)+i                                            
      indic(ij)=ijs                                                     
      jndic(ij)=ijn                                                     
      lndic(ij)=ijm                                                     
50    continue                                                          
51    continue                                                          
c     ij=0                                                              
c                                                                       
700   format(1x,20i6)                                                   
c     write(2,*) 'ndeb',(ndeb(i),i=1,npair)
c     write(2,*) 'indic, jndic sur 2 lignes consecutives'
c     do 60 i=1,norb                                                    
c     ij=num(i)+1                                                       
c     ij2=num(i+1)                                                      
c     write(2,70) (indic(j),j=ij,ij2)                                   
c     write(2,70) (jndic(j),j=ij,ij2)                                   
c     write(2,70) (lndic(j),j=ij,ij2)                                   
c60    write(2,70)                                                      
70    format(1x,30i4)                                                   
      return                                                            
      end                                                               
