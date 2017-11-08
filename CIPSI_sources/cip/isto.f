      subroutine isto(val,vil,nespf,fmul,metat,teseff)                  
c                                                                       
c                                                                       
c        subroutine che fa l'istogramma dei determinanti secondo        
c        il loro contributo alla funzione d'onda                        
c                                                                       
c                                                                       
      implicit real*8(a-h,o-x,z),logical*1(y)                           
      include 'pshf.prm'
      common/hist/nom(10),nam(10,metz),aval(10),tabmp(10,metz),         
     1tau,semp(metz),sevp(metz),tabvp(10,metz),nclass,ymoyen            
      dimension vil(metz)                                               
      if(val.lt.teseff) return                                          
      nc=nclass-1                                                       
      fmul=fmul+1.                                                      
      do 1 i=1,nc                                                       
      if(val.lt.aval(i))go to 10                                        
 1    continue                                                          
      nom(nclass)=nom(nclass)+1+nespf                                   
      do 2 m=1,metat                                                    
      if(vil(m).lt.teseff) go to 2                                      
      tabmp(nclass,m)=tabmp(nclass,m)+semp(m)*fmul                      
      tabvp(nclass,m)=tabvp(nclass,m)+sevp(m)*fmul                      
2     continue                                                          
      go to 30                                                          
 10   nom(i)=nom(i)+1+nespf                                             
      do 3 m=1,metat                                                    
      if(vil(m).lt.teseff) go to 3                                      
      tabmp(i,m)=tabmp(i,m)+semp(m)*fmul                                
      tabvp(i,m)=tabvp(i,m)+sevp(m)*fmul                                
3     continue                                                          
 30   continue                                                          
      do 5 m=1,metat                                                    
      if(vil(m).lt.teseff) go to 5                                      
      do 4 i=1,nc                                                       
      if(vil(m).lt.aval(i)) go to 20                                    
 4    continue                                                          
      nam(nclass,m)=nam(nclass,m)+nespf+1                               
      go to 5                                                           
 20   nam(i,m)=nam(i,m)+nespf+1                                         
 5    continue                                                          
      return                                                            
      end                                                               
