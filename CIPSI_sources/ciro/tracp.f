      function tracp(a,b,n)                                             
      implicit real*8(a-h,o-z)                                          
      include 'pshf.prm'
c      dimension a(1),b(1)                                              
      dimension a(doa**2),b(doas)
      data zero/0.d0/                                                   
      tracp=zero                                                        
      do 10 i=1,n                                                       
      do 20 k=1,n                                                       
      kia=(k-1)*n+i                                                     
      if(k.ge.i) then                                                   
      kib=k*(k-1)/2+i                                                   
      else                                                              
      kib=i*(i-1)/2+k                                                   
      endif                                                             
   20 tracp=tracp+a(kia)*b(kib)                                         
   10 continue                                                          
      return                                                            
      end                                                               
