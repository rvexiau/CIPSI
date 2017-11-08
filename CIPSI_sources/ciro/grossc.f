      subroutine grossc(a,b,n)                                          
      implicit real*8 (a-h,o-z)                                         
      dimension a(2),b(2)                                               
      data zero /0.0d+00/                                               
      ia(i)=i*(i-1)/2                                                   
      do 200 i=1,n                                                      
      dum=zero                                                          
      do 100 j=1,n                                                      
      ij=ia(i)+j                                                        
      if(j.gt.i) ij=ia(j)+i                                             
  100 dum=dum+a(ij)                                                     
  200 b(i)=dum                                                          
      return                                                            
      end                                                               
