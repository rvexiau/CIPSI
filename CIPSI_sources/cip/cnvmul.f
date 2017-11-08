      function cnvmul(fmul)                                             
      implicit real*8(a-h,o-z)                                          
      delta=1.+4.*fmul                                                  
      sol=(-1.+dsqrt(delta))/2.                                         
      cnvmul=sol                                                        
      return                                                            
      end                                                               
