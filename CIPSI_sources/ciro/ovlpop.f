      subroutine ovlpop(a,b,n)                                          
      implicit real*8 (a-h,o-z)                                         
      dimension a(2),b(2)                                               
      do 100 i=1,n                                                      
  100 a(i)=a(i)*b(i)                                                    
      return                                                            
      end                                                               
