      subroutine escriu(av,nf,nc,n)                                     
      implicit real*8 (a-h,o-z)                                         
      dimension av(n,n)                                                 
      character*4 line(31)                                              
      do 888 i=1,31                                                     
888   line(i)='----'                                                    
      ka=1                                                              
      kc=10                                                             
10    kb=min0(kc,nc)                                                    
      write(6,50) (i,i=ka,kb)                                           
      na=3*(kb-ka+1)+1                                                  
      write(6,60) (line(k),k=1,na)                                      
      do 30 i=1,nf                                                      
      write(6,80) i,(av(i,j),j=ka,kb)                                   
30    continue                                                          
      if(kb.eq.nc) return                                               
      ka=kc+1                                                           
      kc=kc+10                                                          
      go to 10                                                          
50    format(/,2x,10(7x,i5))                                            
60    format(1x,31a4)                                                   
80    format(1x,i3,2x,10f12.6)
      end                                                               
