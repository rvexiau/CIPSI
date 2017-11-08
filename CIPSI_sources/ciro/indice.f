       function indice(i,j)                                             
       if(i.gt.j)go to 2                                                
       ij=j*(j-1)/2+i                                                   
       indice=ij                                                        
       return                                                           
2      ij=i*(i-1)/2+j                                                   
       indice=ij                                                        
       return                                                           
       end                                                              
