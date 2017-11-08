      subroutine nomf(fil,prefix,racine)                                
      character*1 prefix(80),racine(10),blanc,fil(90)                   
      blanc= ' '                                                        
      do 1 i=1,80                                                       
	 if(prefix(i).eq.blanc)then                                            
	    do 2 j=1,i-1                                                       
	       fil(j)=prefix(j)                                                
 2          continue                                                    
	    k=i-1                                                              
	    goto 3                                                             
         endif                                                          
 1    continue                                                          
 3    continue                                                          
      do 4 i=1,10                                    
	 if(racine(i).eq.blanc) then                                           
	    do 5 j=1,i-1                                                       
	       fil(k+j)=racine(j)                                              
 5          continue                                                    
	    k=k+i-1                                                            
	    goto 6                                                             
         endif                                                          
 4    continue                                                          
 6    continue      
      return                                                            
      end
