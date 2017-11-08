      subroutine atpop(a,b,nat)                                         
      implicit real*8 (a-h,o-z)                                         
      include 'pshf.prm'
      common/atlim/limlow(dc),limsup(dc)                                
      dimension a(2),b(2)                                               
      data zero /0.0d+00/                                               
      ia(i)=i*(i-1)/2                                                   
      do 200 i=1,nat                                                    
      do 200 j=1,nat                                                    
      i1=limlow(i)                                                      
      i2=limsup(i)                                                      
      j1=limlow(j)                                                      
      j2=limsup(j)                                                      
      dum=zero                                                          
      do 100 k=i1,i2                                                    
      do 100 l=j1,j2                                                    
      kl=ia(k)+l                                                        
      if(l.gt.k) kl=ia(l)+k                                             
  100 dum=dum+a(kl)                                                     
      ij=ia(i)+j                                                        
      if(j.gt.i) ij=ia(j)+i                                             
  200 b(ij)=dum                                                         
      return                                                            
      end                                                               
