      subroutine initab(tabmul)                                         
      implicit real*8(a-h,o-z)                                          
      character*10 tabmul                                               
      dimension tabmul(10)                                              
      tabmul(1)='singulet'                                              
      tabmul(2)='doublet'                                               
      tabmul(3)='triplet'                                               
      tabmul(4)='quadruplet'                                            
      tabmul(5)='quintuplet'                                            
      tabmul(6)='sextuplet'                                             
      tabmul(7)='septuplet'                                             
      tabmul(8)='octuplet'                                              
      tabmul(9)='nonuplet'                                              
      tabmul(10)='---'                                                  
      return                                                            
      end                                                               
