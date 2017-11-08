      subroutine initsm(tabsym,nsym)                                    
      implicit real*8(a-h,o-z)                                          
      common/grou/group                                                 
      character*8 group                                                 
      character*7 tabsym                                                
      dimension tabsym(10)                                              
      if (nsym.eq.5) then                                               
         group='cinfv'                                                  
         tabsym(1)='s '                                                 
         tabsym(2)='px'                                                 
         tabsym(3)='py'                                                 
         tabsym(4)='dx2-y2'                                             
         tabsym(5)='dxy'                                                
      else if (nsym.eq.1) then                                          
         group='c1'                                                     
         tabsym(1)='a'                                                  
      else if (nsym.eq.2) then                                          
         group='cs'                                                     
         tabsym(1)='a'' '                                               
         tabsym(2)='a'''''                                              
      else if (nsym.eq.4) then                                          
         group='c2v'                                                    
         tabsym(1)='a1 '                                                
         tabsym(2)='a2 '                                                
         tabsym(3)='b2 '                                                
         tabsym(4)='b1 '                                                
      else if (nsym.eq.8) then                                          
         group='d2h'                                                    
         tabsym(1)='ag '                                                
         tabsym(2)='au '                                                
         tabsym(3)='b2u'                                                
         tabsym(4)='b2g'                                                
         tabsym(5)='b1g'                                                
         tabsym(6)='b1u'                                                
         tabsym(7)='b3u'                                                
         tabsym(8)='b3g'                                                
      else if (nsym.eq.10) then                                         
         group='dinfh'                                                  
         tabsym(1)='sg'                                                 
         tabsym(2)='su'                                                 
         tabsym(3)='pgx'                                                
         tabsym(4)='pux'                                                
         tabsym(5)='pgy'                                                
         tabsym(6)='puy'                                                
         tabsym(7)='dgx2-y2'                                            
         tabsym(8)='dux2-y2'                                            
         tabsym(9)='dgxy'                                               
         tabsym(10)='duxy'                                              
      else                                                              
         do i=1,10                                                      
            tabsym(i)='...'                                             
         enddo                                                          
      endif                                                             
      return                                                            
      end                                                               
