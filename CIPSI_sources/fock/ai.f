      function ai(i,j,k,l,bijkls,bijkld,ydp)                    
      implicit real*8(a-h,o-x,z),logical*1(y)                           
      include 'pshf.prm'
      real*4 bijkls(*)
      dimension bijkld(*)
c                                                                       
      common/truc/indic(doa*(doa+1)/2),jndic(doa*(doa+1)/2),
     & kt(10000),nbo(99),ndeb(500),nad(10000),
     & lndic(doa*(doa+1)/2)                                   
c     integer*2 indic,jndic,nad,kt,lndic
      integer*2 indic,jndic,nad,kt
      num(i)=i*(i-1)/2
      if(i.ge.j) go to 10                                               
      ij=num(j)+i                                                       
      go to 15                                                          
10    ij=num(i)+j                                                       
15    if(k.ge.l) go to 20                                               
      kl=num(l)+k                                                       
      go to 25                                                          
20    kl=num(k)+l                                                       
25    ijs=indic(ij)                                                     
       kls=indic(kl)                                                    
      if(ijs.le.kls) go to 30                                           
       ijkls=num(ijs)+kls                                               
      go to 35                                                          
30    ijkls=num(kls)+ijs                                                
      nij=ij                                                            
      ij=kl                                                             
      kl=nij                                                            
      ijs=kls                                                           
      kls=indic(kl)                                                     
35    continue                                                          
      if(nad(ijkls)) 40,300,45                                          
40    ysig=.true.                                                       
      nadi=-nad(ijkls)                                                  
      go to 50                                                          
45    ysig=.false.                                                      
      nadi=nad(ijkls)                                                   
50    ktyp=kt(ijkls)                                                    
      nij=jndic(ij)                                                     
      go to (60,70,75,80,85,90,95,100,105,62),ktyp                      
62    nij=lndic(ij)                                                     
      nkl=lndic(kl)                                                     
      go to 63                                                          
60    nij=jndic(ij)                                                     
      nkl=jndic(kl)                                                     
63    if(nij.le.nkl) go to 65                                           
      ijkl=num(nij)+nkl                                                 
      go to 200                                                         
65    ijkl=num(nkl)+nij                                                 
      go to 200                                                         
70    nij=jndic(ij)                                                     
      nkl=jndic(kl)                                                     
      ijkl=nbo(kls)*(nij-1)+nkl                                         
      go to 200                                                         
75    nij=lndic(ij)                                                     
      nkl=jndic(kl)                                                     
      ijkl=nbo(kls)*(nij-1)+nkl                                         
      go to 200                                                         
80    nij=jndic(ij)                                                     
      nkl=lndic(kl)                                                     
      ijkl=nbo(kls)*(nij-1)+nkl                                         
      go to 200                                                         
85    nij=lndic(ij)                                                     
      nkl=lndic(kl)                                                     
      ijkl=nbo(kls)*(nij-1)+nkl                                         
      go to 200                                                         
90    nij=jndic(ij)                                                     
      nkl=jndic(kl)                                                     
      ijkl=nbo(ijs)*(nkl-1)+nij                                         
      go to 200                                                         
95    nij=lndic(ij)                                                     
      nkl=jndic(kl)                                                     
      ijkl=nbo(ijs)*(nkl-1)+nij                                         
      go to 200                                                         
100   nij=jndic(ij)                                                     
      nkl=lndic(kl)                                                     
      ijkl=nbo(ijs)*(nkl-1)+nij                                         
      go to 200                                                         
105   nij=lndic(ij)                                                     
      nkl=lndic(kl)                                                     
      ijkl=nbo(ijs)*(nkl-1)+nij                                         
      go to 200                                                         
200   continue
      if(.not.ydp) ai=bijkls(ndeb(nadi)+ijkl)                           
      if(ydp) ai=bijkld(ndeb(nadi)+ijkl)                                
      if(ysig) ai=-ai                                                   
      return                                                            
300   ai=0.d0                                                           
      return                                                            
      end                                                               
